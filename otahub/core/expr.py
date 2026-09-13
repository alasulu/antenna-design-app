"""Safe evaluation of the formula strings carried in the archetype specs.

Specs are data, not code: every ``expr`` field is a restricted Python
expression evaluated against a whitelisted namespace.  Anything outside the
whitelist is a hard error rather than a silent fallback, because a formula
that quietly evaluates to the wrong thing is worse than one that refuses.
"""
from __future__ import annotations

import ast
import cmath
import math
from typing import Any, Mapping

import numpy as np
from scipy import special as _sp

from .constants import C0, EPS0, ETA0, MU0

#: Euler-Mascheroni constant, which appears throughout wire-antenna theory.
EULER_GAMMA = 0.5772156649015329


def _si(x):
    """Sine integral Si(x) = int_0^x sin(t)/t dt."""
    return _sp.sici(x)[0]


def _ci(x):
    """Cosine integral Ci(x) = -int_x^inf cos(t)/t dt."""
    return _sp.sici(x)[1]

# Node types the evaluator will walk. Anything else raises.
_ALLOWED_NODES: tuple[type[ast.AST], ...] = (
    ast.Expression, ast.Constant, ast.Name, ast.Load,
    ast.BinOp, ast.UnaryOp, ast.Call, ast.Attribute,
    ast.Add, ast.Sub, ast.Mult, ast.Div, ast.FloorDiv, ast.Mod, ast.Pow,
    ast.USub, ast.UAdd,
    ast.Compare, ast.Eq, ast.NotEq, ast.Lt, ast.LtE, ast.Gt, ast.GtE,
    ast.IfExp, ast.BoolOp, ast.And, ast.Or, ast.Not,
    ast.Tuple, ast.List, ast.Subscript, ast.Index if hasattr(ast, "Index") else ast.Load,
)

# numpy attributes a spec may reach through ``np.``
_NP_WHITELIST = frozenset("""
    sqrt exp log log10 log2 sin cos tan arcsin arccos arctan arctan2 sinh cosh tanh
    abs pi e inf deg2rad rad2deg degrees radians power maximum minimum clip
    real imag angle conj sign floor ceil round where array linspace arange
    sum prod mean max min argmax argmin interp trapezoid
""".split())

#: Constants and functions every spec expression may use unqualified.
BASE_NAMESPACE: dict[str, Any] = {
    "c": C0, "c0": C0, "eps0": EPS0, "mu0": MU0, "eta0": ETA0,
    "pi": math.pi, "e": math.e, "j": 1j, "inf": float("inf"),
    "sqrt": np.sqrt, "exp": np.exp, "log": np.log, "ln": np.log,
    "log10": np.log10, "log2": np.log2,
    "sin": np.sin, "cos": np.cos, "tan": np.tan,
    "asin": np.arcsin, "acos": np.arccos, "atan": np.arctan, "atan2": np.arctan2,
    "arcsin": np.arcsin, "arccos": np.arccos, "arctan": np.arctan,
    "sinh": np.sinh, "cosh": np.cosh, "tanh": np.tanh,
    "abs": np.abs, "sign": np.sign, "floor": np.floor, "ceil": np.ceil,
    "deg": np.rad2deg, "rad": np.deg2rad,
    "min": min, "max": max, "round": round,
    # Special functions that closed-form antenna theory genuinely requires:
    # the exact dipole radiation resistance is written in Si/Ci, and circular
    # apertures and loops in Bessel functions. Without these a spec has to
    # fall back on fitted approximations.
    "Si": _si, "Ci": _ci, "gamma_e": EULER_GAMMA, "euler_gamma": EULER_GAMMA,
    "j0": _sp.j0, "j1": _sp.j1, "jv": _sp.jv, "jn": _sp.jv,
    "struve": _sp.struve, "ellipk": _sp.ellipk, "ellipe": _sp.ellipe,
    "db10": lambda x: 10.0 * np.log10(x),
    "db20": lambda x: 20.0 * np.log10(x),
    "undb10": lambda x: 10.0 ** (np.asarray(x) / 10.0),
    "np": np,
}


class ExprError(ValueError):
    """Raised when a spec expression is malformed, unsafe, or unevaluable."""


class _Validator(ast.NodeVisitor):
    """Reject any construct outside the whitelist before evaluation."""

    def __init__(self, allowed_names: frozenset[str]) -> None:
        self.allowed_names = allowed_names
        self.used_names: set[str] = set()

    def generic_visit(self, node: ast.AST) -> None:
        if not isinstance(node, _ALLOWED_NODES):
            raise ExprError(
                f"disallowed syntax {type(node).__name__!r} in spec expression"
            )
        super().generic_visit(node)

    def visit_Name(self, node: ast.Name) -> None:
        if node.id not in self.allowed_names:
            raise ExprError(
                f"unknown symbol {node.id!r}; declare it as a parameter or add it "
                f"to the base namespace"
            )
        self.used_names.add(node.id)

    def visit_Attribute(self, node: ast.Attribute) -> None:
        # Only ``np.<whitelisted>`` is reachable; blocks __globals__ style escapes.
        if not (isinstance(node.value, ast.Name) and node.value.id == "np"):
            raise ExprError("attribute access is only permitted on 'np'")
        if node.attr not in _NP_WHITELIST:
            raise ExprError(f"numpy attribute {node.attr!r} is not whitelisted")
        self.used_names.add("np")

    def visit_Call(self, node: ast.Call) -> None:
        if node.keywords:
            raise ExprError("keyword arguments are not permitted in spec expressions")
        self.visit(node.func)
        for arg in node.args:
            self.visit(arg)


def compile_expr(expr: str, param_names: frozenset[str] = frozenset()) -> ast.Expression:
    """Parse and validate `expr`, returning the AST ready for evaluation."""
    try:
        tree = ast.parse(expr.strip(), mode="eval")
    except SyntaxError as exc:
        raise ExprError(f"cannot parse {expr!r}: {exc}") from exc
    allowed = frozenset(BASE_NAMESPACE) | param_names
    _Validator(allowed).visit(tree)
    return tree


def referenced_symbols(expr: str, param_names: frozenset[str] = frozenset()) -> set[str]:
    """Return the parameter symbols an expression actually reads.

    Used to cross-check a spec's declared ``depends_on`` — a wrong dependency
    list would silently break the synthesis topological sort.
    """
    tree = ast.parse(expr.strip(), mode="eval")
    validator = _Validator(frozenset(BASE_NAMESPACE) | param_names)
    validator.visit(tree)
    return {n for n in validator.used_names if n not in BASE_NAMESPACE}


def evaluate(expr: str, variables: Mapping[str, Any]) -> Any:
    """Evaluate a spec expression against `variables`."""
    tree = compile_expr(expr, frozenset(variables))
    namespace = dict(BASE_NAMESPACE)
    namespace.update(variables)
    code = compile(tree, filename="<spec>", mode="eval")
    try:
        return eval(code, {"__builtins__": {}}, namespace)  # noqa: S307 - validated above
    except ExprError:
        raise
    except ZeroDivisionError as exc:
        raise ExprError(f"division by zero evaluating {expr!r}") from exc
    except Exception as exc:  # noqa: BLE001 - surface the spec's fault with context
        raise ExprError(f"failed evaluating {expr!r}: {type(exc).__name__}: {exc}") from exc

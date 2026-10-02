"""Design values that remember how they were computed.

The CST and HFSS scripts declare every design value as a variable, but they
used to write the solids with the numbers, so editing a variable moved nothing
and the simulator's optimiser could not drive the design. An `Expr` is a float -
every numeric use of a model (the 3-D mesh, the checks, the notes) sees exactly
the number it always did - that also carries its expression in those
variables, built up by the builders' own arithmetic. The renderers write the
expression, in their own units and syntax, where there is one, and the number
where there is not.

An expression is kept as a sum of coefficient-times-atom terms plus a constant,
an atom being a variable, or a product, quotient or function that does not
reduce further: a symmetric brick's size comes out as W, not W/2 - (-W/2).
Every value carries its physical dimension (powers of length and time), so a
constant can be written in the renderer's own units - 1.5 mm is "1.5" to CST,
which works in millimetres, and "1.5mm" to HFSS, which checks units.
"""
from __future__ import annotations

import math
from decimal import Decimal
from typing import Callable, NamedTuple

#: dimensions, as (length, time) exponents
NONE, LENGTH = (0, 0), (1, 0)


def unit_dimension(unit: str) -> tuple[int, int]:
    """The dimension of a variable the spec declares in `unit`. Only metres and
    hertz are declared as lengths and frequencies; angles and everything else
    are pure numbers."""
    return {"m": (1, 0), "Hz": (0, -1)}.get((unit or "").strip(), NONE)


class Expr(float):
    """A float that also knows its expression in the design's variables.

    `lin` is (terms, constant, dimension), the terms a tuple of (atom,
    coefficient). Arithmetic with plain numbers and other Exprs gives an Expr
    whose float value is exactly what the plain arithmetic gives.
    """
    __slots__ = ("lin",)

    def __new__(cls, value: float, lin: tuple):
        self = super().__new__(cls, value)
        self.lin = lin
        return self

    def __reduce__(self):                       # copy and pickle keep the expression
        return Expr, (float(self), self.lin)

    @property
    def dim(self) -> tuple[int, int]:
        return self.lin[2]

    def __add__(self, other):
        return _add(self, other, float(self) + float(other), 1) if _number(other) else NotImplemented

    def __radd__(self, other):
        return _add(other, self, float(other) + float(self), 1) if _number(other) else NotImplemented

    def __sub__(self, other):
        return _add(self, other, float(self) - float(other), -1) if _number(other) else NotImplemented

    def __rsub__(self, other):
        return _add(other, self, float(other) - float(self), -1) if _number(other) else NotImplemented

    def __mul__(self, other):
        return _mul(self, other, float(self) * float(other)) if _number(other) else NotImplemented

    def __rmul__(self, other):
        return _mul(other, self, float(other) * float(self)) if _number(other) else NotImplemented

    def __truediv__(self, other):
        return _div(self, other, float(self) / float(other)) if _number(other) else NotImplemented

    def __rtruediv__(self, other):
        return _div(other, self, float(other) / float(self)) if _number(other) else NotImplemented

    def __neg__(self):
        return _result(-float(self), *_scaled(self.lin, -1.0), self.dim)

    def __pos__(self):
        return self

    def __abs__(self):
        return _call("abs", abs, self)


def variable(name: str, value: float, unit: str) -> Expr:
    """The design variable `name`, whose value is `value`, declared in `unit`."""
    return Expr(value, (((("var", name, unit), 1.0),), 0.0, unit_dimension(unit)))


def constant(value: float, dim: tuple[int, int]) -> Expr:
    """A physical constant with a dimension, such as the speed of light. A bare
    number in a builder's arithmetic counts as a pure factor (or, added to a
    quantity, as that quantity's kind)."""
    return Expr(value, ((), value, dim))


def sqrt(x):
    """math.sqrt, followed by the scripts when what is under the root has a whole
    dimension's square root (a pure number, an area)."""
    value = math.sqrt(x)
    if not isinstance(x, Expr) or x.dim[0] % 2 or x.dim[1] % 2:
        return value
    return _result(value, ((("call", "sqrt", x.lin), 1.0),), 0.0, (x.dim[0] // 2, x.dim[1] // 2))


def sin(x):
    return _call("sin", math.sin, x)


def cos(x):
    return _call("cos", math.cos, x)


def maximum(*values):
    """max(), in a form the scripts can follow: max(a, b) = (a + b + |a - b|)/2,
    which every expression language can write, so a feed gap set by the larger
    of two rules switches rule in the simulator as it does here. The value is
    max()'s own."""
    return _extreme(max, 1, values)


def minimum(*values):
    """min(), as (a + b - |a - b|)/2; see maximum."""
    return _extreme(min, -1, values)


# ------------------------------------------------------------------ rendering

class Dialect(NamedTuple):
    """How one simulator writes an expression."""
    #: (name, unit) -> the token for that variable, holding the value the builder saw
    var: Callable[[str, str], str]
    #: (SI value, dimension) -> a constant in the simulator's units
    num: Callable[[float, tuple], str]
    #: sqrt, sin, cos, abs -> the simulator's function names
    calls: dict


def expression(value, dialect: Dialect, dim: tuple[int, int] = LENGTH) -> str | None:
    """`value` as an expression in `dialect`, or None when it carries none (or
    is not of dimension `dim`) and the caller should write the number."""
    if not isinstance(value, Expr) or value.dim != dim:
        return None
    return _lin_text(value.lin, dialect)[0]


def number(value: float) -> str:
    """A number to 15 significant figures - closer than any simulator resolves
    - never in exponent form, which a unit suffix would not survive ("1e-05mm"),
    and without a trailing ".0"."""
    if value == 0:
        return "0"
    out = format(Decimal(f"{value:.15g}"), "f")
    return out.rstrip("0").rstrip(".") if "." in out else out


# ------------------------------------------------------------------ internals

_ONE = ((), 1.0, NONE)                          # the numerator of a reciprocal


def _number(x) -> bool:
    return isinstance(x, (int, float))


def _plus(a, b):
    return a[0] + b[0], a[1] + b[1]


def _minus(a, b):
    return a[0] - b[0], a[1] - b[1]


def _lin(x, dim):
    """x's linear form; a bare number is a constant of dimension `dim`."""
    return x.lin if isinstance(x, Expr) else ((), float(x), dim)


def _result(value, terms, const, dim):
    """An Expr, unless nothing variable is left: then the number itself."""
    return Expr(value, (terms, const, dim)) if terms else value


def _scaled(lin, k, divide=False):
    """(terms, constant) of lin times k, or divided by k."""
    op = (lambda c: c / k) if divide else (lambda c: c * k)
    return tuple((atom, op(c)) for atom, c in lin[0] if op(c) != 0.0), op(lin[1])


def _add(a, b, value, sign):
    dim = a.dim if isinstance(a, Expr) else b.dim
    la, lb = _lin(a, dim), _lin(b, dim)
    if la[2] != lb[2]:                          # not one kind of quantity: keep the number
        return value
    merged = dict(la[0])
    for atom, c in lb[0]:
        merged[atom] = merged.get(atom, 0.0) + sign * c
    terms = tuple((atom, c) for atom, c in merged.items() if c != 0.0)
    return _result(value, terms, la[1] + sign * lb[1], dim)


def _mul(a, b, value):
    la, lb = _lin(a, NONE), _lin(b, NONE)
    dim = _plus(la[2], lb[2])
    if not la[0]:
        return _result(value, *_scaled(lb, la[1]), dim)
    if not lb[0]:
        return _result(value, *_scaled(la, lb[1]), dim)
    return _result(value, ((("mul", la, lb), 1.0),), 0.0, dim)


def _div(a, b, value):
    la, lb = _lin(a, NONE), _lin(b, NONE)
    dim = _minus(la[2], lb[2])
    if not lb[0]:
        return _result(value, *_scaled(la, lb[1], divide=True), dim)
    if not la[0]:                               # k / b: k times a reciprocal
        reciprocal = (((("div", _ONE, lb), 1.0),), 0.0, dim)
        return _result(value, *_scaled(reciprocal, la[1]), dim)
    return _result(value, ((("div", la, lb), 1.0),), 0.0, dim)


def _call(name, fn, x):
    value = fn(float(x))
    if not isinstance(x, Expr) or (name != "abs" and x.dim != NONE):
        return value
    return _result(value, ((("call", name, x.lin), 1.0),), 0.0, x.dim if name == "abs" else NONE)


def _extreme(pick, sign, values):
    best = values[0]
    for x in values[1:]:
        value = pick(float(best), float(x))
        if isinstance(best, Expr) or isinstance(x, Expr):
            form = (best + x + sign * abs(best - x)) / 2.0
            best = Expr(value, form.lin) if isinstance(form, Expr) else value
        else:
            best = value
    return best


def _atom_dim(atom):
    kind = atom[0]
    if kind == "var":
        return unit_dimension(atom[2])
    if kind in ("mul", "div"):
        return (_plus if kind == "mul" else _minus)(atom[1][2], atom[2][2])
    inner = atom[2][2]                          # a call
    return {"sqrt": (inner[0] // 2, inner[1] // 2), "abs": inner}.get(atom[1], NONE)


# precedence: 1 a sum or a leading minus, 2 a product or quotient, 3 atomic

def _wrap(s, prec, need):
    return f"({s})" if prec < need else s


def _prec(s):
    """The precedence of a token or constant from the dialect ("1.5mm" is one)."""
    if s.startswith("-") or any(ch in s[1:] for ch in "+-"):
        return 1
    return 2 if any(ch in s for ch in "*/") else 3


def _lin_text(lin, d):
    terms, const, dim = lin
    parts = [_term_text(atom, c, dim, d) for atom, c in terms]
    if const or not parts:
        s = d.num(const, dim)
        parts.append((s, _prec(s)))
    if len(parts) == 1:
        return parts[0]
    out = parts[0][0]
    for s, _ in parts[1:]:
        out += s if s.startswith("-") else "+" + s
    return out, 1


def _term_text(atom, c, dim, d):
    cdim = _minus(dim, _atom_dim(atom))
    if atom[0] == "div" and atom[1] == _ONE:    # c / b
        k, (den, p) = d.num(c, cdim), _lin_text(atom[2], d)
        return f"{k}/{_wrap(den, p, 3)}", (1 if k.startswith("-") else 2)
    body, p = _atom_text(atom, d)
    if cdim == NONE:
        if c == 1.0:
            return body, p
        if c == -1.0:
            return "-" + _wrap(body, p, 2), 1
        n = 1.0 / c                             # W/2 reads better than 0.5*W
        if abs(n) >= 2 and abs(n - round(n)) <= 1e-12 * abs(n):
            return ("-" if n < 0 else "") + f"{_wrap(body, p, 2)}/{abs(round(n))}", (1 if n < 0 else 2)
    k = d.num(c, cdim)
    return f"{k}*{_wrap(body, p, 2)}", (1 if k.startswith("-") else 2)


def _atom_text(atom, d):
    kind = atom[0]
    if kind == "var":
        token = d.var(atom[1], atom[2])
        return token, _prec(token)
    if kind == "call":
        return f"{d.calls[atom[1]]}({_lin_text(atom[2], d)[0]})", 3
    (left, lp), (right, rp) = _lin_text(atom[1], d), _lin_text(atom[2], d)
    if kind == "mul":
        return f"{_wrap(left, lp, 2)}*{_wrap(right, rp, 2)}", 2
    return f"{_wrap(left, lp, 2)}/{_wrap(right, rp, 3)}", 2

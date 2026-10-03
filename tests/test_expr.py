"""The expression evaluator is the trust boundary: specs are data, not code."""
import math

import pytest

from otahub.core.expr import ExprError, evaluate, referenced_symbols


def test_basic_arithmetic_and_constants():
    assert evaluate("c", {}) == pytest.approx(2.99792458e8)
    assert evaluate("eta0", {}) == pytest.approx(376.730313412, rel=1e-9)
    assert evaluate("0.4788 * c / f0", {"f0": 300e6}) == pytest.approx(0.478469, rel=1e-4)


def test_complex_impedance_round_trips():
    z = evaluate("73.1 + 42.5j", {})
    assert isinstance(z, complex)
    assert z.real == pytest.approx(73.1) and z.imag == pytest.approx(42.5)


def test_numpy_whitelist_allows_sqrt_but_not_arbitrary_attributes():
    assert evaluate("np.sqrt(9)", {}) == pytest.approx(3.0)
    with pytest.raises(ExprError, match="not whitelisted"):
        evaluate("np.load('x')", {})


@pytest.mark.parametrize("hostile", [
    "__import__('os').system('id')",
    "(1).__class__.__bases__",
    "open('/etc/passwd').read()",
    "[x for x in range(3)]",
    "lambda: 1",
])
def test_hostile_expressions_are_rejected(hostile):
    with pytest.raises(ExprError):
        evaluate(hostile, {})


def test_unknown_symbol_is_an_error_not_a_zero():
    # Silently treating an unknown symbol as 0 would produce a plausible,
    # wrong antenna. It must fail loudly instead.
    with pytest.raises(ExprError, match="unknown symbol"):
        evaluate("length_that_was_never_defined * 2", {})


def test_referenced_symbols_reports_actual_dependencies():
    deps = referenced_symbols("W / h + np.sqrt(eps_r)", frozenset({"W", "h", "eps_r"}))
    assert deps == {"W", "h", "eps_r"}


def test_division_by_zero_is_reported_with_context():
    with pytest.raises(ExprError, match="division by zero"):
        evaluate("c / f0", {"f0": 0})


def test_a_domain_error_is_reported_in_numpys_words():
    """numpy's warning for log10(-5) cannot run in the evaluator's builtin-free frame,
    and surfaced as "KeyError: '__import__'" - a message about nothing."""
    from otahub.core.expr import ExprError, evaluate
    for x, says in ((-5.0, "invalid value"), (0.0, "divide by zero")):
        with pytest.raises(ExprError, match=says):
            evaluate("log10(x)", {"x": x})

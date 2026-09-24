"""Top-loaded monopole: the top current derived from the hat, at last.

Since session 5 this spec carried its hat radius and its top-current ratio as
two independent inputs, and said plainly why: relating them "needs a numerical
solve this spec does not attempt". Junctions made that solve possible - a
vertical with a hat of radial wires, joined where they meet, over infinite
ground by image theory - and beta is now derived from the hat. Supplying it
still overrides, so every design that set it explicitly is unchanged.

What the solve also settled: the spec's linear-taper radiation resistance was
never the weak link. Given the right beta it is within 2.4% of the radiated
power for every hat tried. It only ever needed the right beta.
"""
from __future__ import annotations

import functools
import math

import numpy as np
import pytest

from otahub.num import mom

C0 = 2.99792458e8


@functools.lru_cache(maxsize=32)
def _solve(h, R, N, aw):
    m, feed, top = mom.top_hat_monopole(h, R, N, aw)
    sol = mom.solve(m, feed)
    nv = 39
    z = np.array([m.node_of(k)[2] for k in range(nv)])
    sel = (z > 0.25 * h) & (z < 0.85 * h)
    c = np.polyfit(z[sel], np.abs(sol.currents[:nv][sel]), 1)
    i0 = np.polyval(c, 0.0)
    return dict(beta=np.polyval(c, h) / i0,
                r=2 * (0.5 * mom.radiated_power(sol, 40, 40)) / i0 ** 2,
                x=sol.input_impedance.imag / 2,
                d=2 * mom.directivity(sol, 40, 40))


def _spec(registry, h=0.05, R=0.01, N=8, aw_l=3.3e-5, **kw):
    lam = C0 / 1e7
    return registry["top_loaded_monopole"].synthesize(
        f0=1e7, h_over_lambda=h, hat_radius_over_lambda=R, hat_radials=N,
        aw=aw_l * lam, **kw)


@pytest.mark.parametrize("h,R,N", [(0.05, 0.01, 8), (0.05, 0.02, 16), (0.03, 0.02, 4),
                                   (0.075, 0.04, 8)])
def test_derived_beta_matches_the_solver(h, R, N, registry):
    want = _solve(h, R, N, 3.3e-5)["beta"]
    got = _spec(registry, h, R, N).get("beta_top")
    assert got == pytest.approx(want, abs=0.03)


@pytest.mark.parametrize("h,R,N", [(0.05, 0.01, 4), (0.05, 0.04, 16), (0.03, 0.02, 8)])
def test_linear_taper_resistance_holds_given_the_right_beta(h, R, N):
    """The model was never the problem; the input was."""
    s = _solve(h, R, N, 3.3e-5)
    rr = 160 * math.pi ** 2 * (h * (1 + s["beta"]) / 2) ** 2
    assert s["r"] == pytest.approx(rr, rel=0.03)


def test_an_explicit_beta_still_overrides(registry):
    d = _spec(registry, beta_top=0.6)
    assert d.get("beta_top") == 0.6
    assert d.metrics["radiation_resistance_ohm"] == pytest.approx(
        160 * math.pi ** 2 * (0.05 * 1.6 / 2) ** 2, rel=1e-9)


def test_the_old_default_beta_did_not_match_the_old_default_hat(registry):
    """0.6 sat beside a 0.01 lambda hat; that hat, with eight radials, gives
    about 0.52, and the resistance the old pairing implied was 11% high."""
    d = _spec(registry)
    assert 0.45 < d.get("beta_top") < 0.58
    old = 160 * math.pi ** 2 * (0.05 * 1.6 / 2) ** 2
    assert old / d.metrics["radiation_resistance_ohm"] > 1.08


def test_beta_moves_the_right_way_with_every_hat_parameter(registry):
    base = _spec(registry).get("beta_top")
    assert _spec(registry, R=0.02).get("beta_top") > base, "bigger hat"
    assert _spec(registry, N=16).get("beta_top") > base, "more radials"
    assert _spec(registry, h=0.075).get("beta_top") < base, "taller whip, same hat"


def test_more_radials_keep_helping_past_the_fitted_range():
    """So a solid disc beats any radial count the fit covers."""
    betas = [_solve(0.05, 0.01, n, 1e-5)["beta"] for n in (8, 16, 32)]
    assert betas[0] < betas[1] < betas[2]
    assert betas[2] - betas[1] > 0.03


def test_directivity_is_still_three(registry):
    for R in (0.01, 0.04):
        assert _solve(0.05, R, 8, 3.3e-5)["d"] == pytest.approx(3.0, abs=0.03)


def test_the_reactance_table_includes_self_resonant_hats(registry):
    """The reactance ships as solved values, not a fit - and includes the hats
    big enough to resonate the whip by themselves, which the beta fit excludes."""
    tab = registry["top_loaded_monopole"].spec.tables["mom_top_hat"]
    xi = tab["columns"].index("input_reactance_ohm")
    xs = [row[xi] for row in tab["rows"]]
    assert min(xs) < -1000 and max(xs) > 0

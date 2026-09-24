"""The slot's resonance, arbitrated by two solvers that share nothing.

By Babinet a slot in a screen resonates where its complementary dipole does,
and that dipole is a wire of radius w/4. `half_wave_slot` used to carry a thin
dipole's 0.4785 lambda and a flat 67 ohm whatever its width; it now takes both
from `resonant_dipole` at the complementary radius. At the default width that
radius is 0.0058 lambda, just past where the dipole fits were made, so the fits
alone are not enough - they are checked here against the exact-kernel method of
moments and against Hallen's equation, each at meshes it can be trusted on.

Hallen uses the reduced kernel, so on the fat wire it is kept to segments of
two radii or more; its own docstring says why.
"""
from __future__ import annotations

import functools
import warnings

import pytest
from scipy.optimize import brentq

from otahub.num import mom
from otahub.num.hallen import hallen_dipole

pytestmark = pytest.mark.slow   # full-wave solves; skip with -m 'not slow'

C0 = 2.99792458e8
F = 3e8
ETA0 = 376.730313412


def _resonance(zfn):
    length = brentq(lambda L: zfn(L).imag, 0.43, 0.50, xtol=1e-6)
    return length, zfn(length).real


@functools.lru_cache(maxsize=8)
def _arbiters(a):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        hallen = [_resonance(lambda L: hallen_dipole(L, a, n=n)) for n in (24, 40)]
        # 42 segments and up: the exact-kernel resonance has converged to 0.05%
        # there, where 22 segments is still 0.1% out
        emfie = [_resonance(lambda L: mom.input_impedance(mom.dipole(L, a, n), exact=True))
                 for n in (41, 81)]
    return hallen + emfie


@pytest.mark.parametrize("w_over_L", [0.02, 0.05])
def test_both_solvers_put_the_slot_where_the_spec_does(registry, w_over_L):
    d = registry["half_wave_slot"].synthesize(f0=F, w_over_L=w_over_L)
    lam = C0 / F
    got = d.parameters["L"] / lam
    for length, _ in _arbiters(round(d.parameters["w"] / 4 / lam, 6)):
        assert got == pytest.approx(length, rel=0.005)


@pytest.mark.parametrize("w_over_L", [0.02, 0.05])
def test_the_complementary_resistance_is_good_to_three_percent(registry, w_over_L):
    """The fit reads 0.6-2.5% HIGH against both solvers (delta-gap resistance
    creeps with the mesh, so neither solver is exact either); what matters is
    that it is nowhere near the 67 ohm the slot used to assume."""
    d = registry["half_wave_slot"].synthesize(f0=F, w_over_L=w_over_L)
    lam = C0 / F
    rd = d.get("Rd_res")
    for _, r in _arbiters(round(d.parameters["w"] / 4 / lam, 6)):
        assert rd == pytest.approx(r, rel=0.03)
        assert r > 72.0


def test_the_old_flat_constants_were_thirteen_percent_high(registry):
    """529.6 ohm from a flat 67; the solvers put the complement at 74.0-75.4 ohm
    at the default width, which is 471-479 ohm for the slot."""
    d = registry["half_wave_slot"].synthesize(f0=F)
    lam = C0 / F
    rs = [r for _, r in _arbiters(round(d.parameters["w"] / 4 / lam, 6))]
    assert ETA0 ** 2 / (4 * max(rs)) > 465.0 and ETA0 ** 2 / (4 * min(rs)) < 485.0
    assert 529.57 / d.metrics["resonant_resistance_ohm"] > 1.10

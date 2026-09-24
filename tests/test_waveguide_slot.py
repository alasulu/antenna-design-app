"""Stevenson's slot conductance, settled from first principles.

The handover had carried "Stevenson's g1 should be cross-checked against the
source" for sessions, and the spec's own note said the wavelength ratio "is
inverted relative to some printings". It was inverted in the spec: the
expression used lambda/lambda_g while the note quoted lambda_g/lambda. Rather
than pick a printing, `otahub.num.waveguide_slot` derives the conductance by
quadrature - reciprocity for the TE10 the slot excites, the slot's one-sided
radiation integrated over a half space, power balance for a shunt element - and
the lambda_g/lambda form is the one it reproduces.

The consequence was not cosmetic. The resonant array sizes its offsets from
g1; with g1 43% low at 10 GHz its twelve WR-90 slots summed to a conductance
of 1.76 instead of 1, a VSWR of 1.76 at a feed designed to be matched.
"""
from __future__ import annotations

import math

import pytest

from otahub.core.constants import ETA0
from otahub.num import waveguide_slot as ws

WR90 = dict(a_wg=0.02286, b_wg=0.01016)
WR28 = dict(a_wg=0.007112, b_wg=0.003556)
BAND = [(WR90, f) for f in (8.2e9, 9.375e9, 10e9, 11e9, 12.4e9)] + \
       [(WR28, f) for f in (26.5e9, 30e9, 35e9, 40e9)]


@pytest.fixture(scope="module")
def slot(registry):
    return registry["waveguide_longitudinal_slot"]


@pytest.mark.parametrize("guide,f0", BAND)
@pytest.mark.parametrize("frac", [0.05, 0.15, 0.25])
def test_the_spec_is_the_first_principles_conductance(slot, guide, f0, frac):
    x1 = frac * guide["a_wg"]
    d = slot.synthesize(f0=f0, x1=x1, **guide)
    want = ws.shunt_conductance(f0, guide["a_wg"], guide["b_wg"], x1)
    assert d.metrics["normalised_conductance"] == pytest.approx(want, rel=2e-4)


@pytest.mark.parametrize("f0,lo,hi", [(8.2e9, 0.60, 0.68), (10e9, 0.40, 0.46), (12.4e9, 0.25, 0.31)])
def test_the_inverted_printing_was_far_too_low(f0, lo, hi):
    """The old expression, lambda/lambda_g, against the derived conductance."""
    a, b = WR90["a_wg"], WR90["b_wg"]
    root = math.sqrt(1 - (2.99792458e8 / (2 * a * f0)) ** 2)            # lambda/lambda_g
    old = 2.09 * root * (a / b) * math.cos(math.pi * root / 2) ** 2
    shortfall = 1 - old / ws.shunt_conductance(f0, a, b, a / 2)
    assert lo < shortfall < hi


def test_coupling_rises_toward_cutoff(slot):
    """The guide's wave impedance climbs toward cutoff and the slot couples
    harder; the inverted form fell toward zero there instead."""
    g = [slot.synthesize(f0=f, x1=0.003, **WR90).metrics["g1_normalised"]
         for f in (12.4e9, 11e9, 10e9, 9e9, 8.2e9)]
    assert all(x < y for x, y in zip(g, g[1:]))


def test_the_twelve_slot_array_is_matched_now(registry):
    """N slots at the designed offset must sum to unity conductance - judged by
    the first-principles model, not by the spec's own g1."""
    d = registry["waveguide_slot_array_resonant"].synthesize(f0=10e9, N=12, **WR90)
    total = 12 * ws.shunt_conductance(10e9, WR90["a_wg"], WR90["b_wg"], d.get("offset"))
    assert total == pytest.approx(1.0, rel=2e-4)
    old_total = 12 * ws.shunt_conductance(10e9, WR90["a_wg"], WR90["b_wg"], 0.0030624)
    gamma = (old_total - 1) / (old_total + 1)
    assert (1 + gamma) / (1 - gamma) > 1.7


def test_the_one_sided_slot_conductance_is_the_induced_emf_dipole_through_babinet():
    """G_ext = 2 R/eta^2 with R = 73.08 ohm: half the two-sided Booker value."""
    assert ws.one_sided_conductance(10e9) == pytest.approx(2 * 73.0796 / ETA0 ** 2, rel=1e-4)


def test_a_real_slot_width_barely_matters():
    a, b = WR90["a_wg"], WR90["b_wg"]
    thin = ws.shunt_conductance(10e9, a, b, 0.25 * a)
    wide = ws.shunt_conductance(10e9, a, b, 0.25 * a, w=1.5e-3)
    assert 0.0 < 1 - wide / thin < 0.005


def test_a_centreline_slot_does_not_couple():
    assert ws.shunt_conductance(10e9, WR90["a_wg"], WR90["b_wg"], 0.0) == pytest.approx(0.0, abs=1e-12)

"""Headless GUI tests.

Run against Qt's offscreen platform so they work in CI and over SSH. These
cover the parts that silently rot - the catalogue filter, the generated form,
and whether every archetype can actually be driven through the UI - rather
than pixel appearance.
"""
from __future__ import annotations

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6", reason="GUI extra not installed")

from PySide6.QtWidgets import QApplication  # noqa: E402

from otahub.core.registry import Registry  # noqa: E402
from otahub.gui.models import (CatalogueFilter, build_catalogue_model,  # noqa: E402
                               default_for, default_frequency, requirement_fields)


@pytest.fixture(scope="session")
def qapp():
    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture(scope="session")
def window(qapp, registry):
    from otahub.gui.app import MainWindow
    return MainWindow(registry)


# ------------------------------------------------------------------ models

def test_catalogue_model_has_one_row_per_family(qapp, registry):
    model = build_catalogue_model(registry)
    assert model.rowCount() == len(registry.families)
    total = sum(model.item(i).rowCount() for i in range(model.rowCount()))
    assert total == len(registry)


def test_filter_narrows_to_matching_family(qapp, registry):
    model = build_catalogue_model(registry)
    proxy = CatalogueFilter(registry)
    proxy.setSourceModel(model)
    proxy.set_needle("horn")
    assert proxy.rowCount() == 1
    assert "horn" in proxy.index(0, 0).data().lower()


def test_filter_restores_everything_when_cleared(qapp, registry):
    model = build_catalogue_model(registry)
    proxy = CatalogueFilter(registry)
    proxy.setSourceModel(model)
    proxy.set_needle("horn")
    proxy.set_needle("")
    assert proxy.rowCount() == len(registry.families)


def test_filter_with_no_matches_hides_every_family(qapp, registry):
    """A family row must not survive on its own - the bug that let the base
    class accept every parent regardless of its children."""
    model = build_catalogue_model(registry)
    proxy = CatalogueFilter(registry)
    proxy.setSourceModel(model)
    proxy.set_needle("no-such-antenna-xyz")
    assert proxy.rowCount() == 0


def test_requirement_fields_put_frequency_first(registry):
    for key in ("rectangular_patch", "half_wave_dipole", "pyramidal_horn"):
        fields = requirement_fields(registry[key])
        assert fields[0].symbol == "f0"


def test_requirement_fields_exclude_derived_and_geometry_outputs(registry):
    roles = {p.role for p in requirement_fields(registry["rectangular_patch"])}
    assert roles <= {"requirement", "assumption", "material"}


def test_default_frequency_lands_inside_the_validity_band(registry):
    for archetype in registry:
        lo, hi = archetype.spec.freq_range_hz
        assert lo <= default_frequency(archetype) <= hi


def test_frequency_field_is_prefilled_even_without_a_declared_typical(registry):
    archetype = registry["half_wave_dipole"]
    f0_param = next(p for p in archetype.spec.parameters if p.symbol == "f0")
    assert f0_param.typical == ""            # the spec declares none
    assert float(default_for(f0_param, archetype)) > 0


# ------------------------------------------------------------------ window

def test_window_builds_with_three_tabs(window):
    assert [window.tabs.tabText(i) for i in range(window.tabs.count())] == \
           ["Catalogue", "Arrays", "Waveguides"]


def test_status_bar_reports_a_clean_catalogue(window, registry):
    assert str(len(registry)) in window.statusBar().currentMessage()


def test_every_archetype_can_be_driven_through_the_form(window, registry):
    """The form is generated from each spec, so this is the test that catches a
    spec whose declared parameters cannot actually produce a design."""
    tab = window.catalogue
    empty = []
    for archetype in registry:
        tab.select_key(archetype.key)
        tab._synthesise()
        if tab.geometry_table.rowCount() + tab.metrics_table.rowCount() == 0:
            empty.append(archetype.key)
    assert not empty, f"produced nothing from form defaults: {empty}"


def test_selecting_an_archetype_populates_its_form(window):
    tab = window.catalogue
    tab.select_key("rectangular_patch")
    assert set(tab._fields) >= {"f0", "eps_r", "h"}
    assert "Rectangular" in tab.title.text()


def test_blank_fields_are_omitted_rather_than_treated_as_zero(window):
    """A blank frequency must not become f0 = 0, which would divide by zero
    deep inside a formula instead of reporting a missing requirement."""
    tab = window.catalogue
    tab.select_key("half_wave_dipole")
    tab._fields["f0"].setText("")
    assert "f0" not in tab.collect()


def test_non_numeric_input_is_rejected_not_coerced(window):
    tab = window.catalogue
    tab.select_key("half_wave_dipole")
    tab._fields["f0"].setText("not a number")
    with pytest.raises(ValueError, match="not a number"):
        tab.collect()


def test_low_confidence_archetypes_are_flagged_in_the_ui(window, registry):
    low = [a for a in registry if a.spec.confidence == "low"]
    assert low, "expected at least one low-confidence archetype to check"
    window.catalogue.select_key(low[0].key)
    assert "low confidence" in window.catalogue.band.text().lower()


def test_pattern_tab_refuses_to_invent_a_pattern(window, registry):
    """Archetypes without a first-principles pattern must say so, not draw a
    plausible-looking guess."""
    from otahub.gui.app import PATTERN_SOURCES
    unsupported = next(a.key for a in registry if a.key not in PATTERN_SOURCES)
    window.catalogue.select_key(unsupported)
    window.catalogue._update_pattern()
    texts = [t.get_text() for t in window.catalogue.pattern_canvas.axes.texts]
    assert any("No closed-form pattern" in t for t in texts)


def test_supported_pattern_actually_plots(window):
    window.catalogue.select_key("half_wave_dipole")
    window.catalogue._update_pattern()
    assert window.catalogue.pattern_canvas.axes.lines


# ------------------------------------------------------------------- tools

def test_array_tab_reports_the_design_sidelobe_level(window):
    tab = window.arrays
    tab.taper.setCurrentText("chebyshev")
    tab.sll.setValue(-25.0)
    tab.n.setValue(16)
    tab.refresh()
    rows = {tab.summary_table.item(r, 0).text(): tab.summary_table.item(r, 1).text()
            for r in range(tab.summary_table.rowCount())}
    assert float(rows["first sidelobe"].split()[0]) == pytest.approx(-25.0, abs=0.1)


def test_array_tab_warns_about_grating_lobes(window):
    tab = window.arrays
    tab.spacing.setValue(1.5)
    tab.scan.setValue(45.0)
    tab.refresh()
    rows = {tab.summary_table.item(r, 0).text(): tab.summary_table.item(r, 1).text()
            for r in range(tab.summary_table.rowCount())}
    assert rows["grating lobe present"] == "YES"
    tab.spacing.setValue(0.5)
    tab.refresh()


def test_waveguide_tab_matches_the_library(window):
    from otahub.waveguides.rectangular import standard
    tab = window.waveguides
    tab.guide.setCurrentText("WR-90")
    tab.refresh()
    rows = {tab.props.item(r, 0).text(): tab.props.item(r, 1).text()
            for r in range(tab.props.rowCount())}
    assert "6.557" in rows["TE10 cutoff"]
    assert tab.modes.rowCount() >= 2
    assert tab.modes.item(0, 0).text() == "TE10"


def test_waveguide_tab_reports_evanescence_below_cutoff(window):
    tab = window.waveguides
    tab.guide.setCurrentText("WR-90")
    tab.freq.setValue(3.0)          # GHz, well below the 6.557 GHz cutoff
    tab.refresh()
    rows = {tab.props.item(r, 0).text(): tab.props.item(r, 1).text()
            for r in range(tab.props.rowCount())}
    assert "CUTOFF" in rows["at this frequency"]
    tab.freq.setValue(10.0)
    tab.refresh()

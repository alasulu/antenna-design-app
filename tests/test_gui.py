"""Headless GUI tests.

Run against Qt's offscreen platform so they work in CI and over SSH. These
cover the parts that silently rot - the catalogue filter, the generated form,
and whether every archetype can actually be driven through the UI - rather
than pixel appearance.
"""
from __future__ import annotations

import math
import os
import re

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

def test_window_builds_with_the_expected_tabs(window):
    assert [window.tabs.tabText(i) for i in range(window.tabs.count())] == \
           ["Catalogue", "Linear arrays", "Planar arrays", "Waveguides"]


def test_status_bar_reports_a_clean_catalogue(window, registry):
    assert str(len(registry)) in window.statusBar().currentMessage()


def test_every_archetype_can_be_driven_through_the_form(window, registry):
    """The form is generated from each spec, so this is the test that catches a
    spec whose declared parameters cannot actually produce a design.

    GEOMETRY specifically, not geometry-or-metrics. Summing the two let a whole
    family pass on nothing but its constant metrics - bandwidth_ratio, a fixed
    directivity - while every dimension came out blank, which is precisely the
    "looks broken" failure the form defaults exist to prevent.
    """
    tab = window.catalogue
    empty = []
    for archetype in registry:
        if not archetype.spec.synthesis:
            continue          # nothing to synthesise: the spec declares no rules
        tab.select_key(archetype.key)
        tab._synthesise()
        if tab.geometry_table.rowCount() == 0:
            empty.append(archetype.key)
    assert not empty, f"produced no geometry from form defaults: {empty}"


def test_a_frequency_requirement_is_prefilled_whatever_it_is_called(window, registry):
    """Keying the prefill to the name "f0" missed every archetype whose
    requirement is `f_low`, so the whole wideband family opened blank."""
    from otahub.gui.models import default_for

    checked = 0
    for archetype in registry:
        for param in archetype.spec.parameters:
            if param.unit == "Hz" and param.role == "requirement":
                assert default_for(param, archetype), (
                    f"{archetype.key}.{param.symbol} opens blank")
                checked += 1
    assert checked > 60, f"only {checked} frequency requirements found"


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


# ------------------------------------------------------------ planar arrays

def _rows(table) -> dict:
    return {table.item(i, 0).text(): table.item(i, 1).text()
            for i in range(table.rowCount())}


def test_planar_tab_is_present(window):
    labels = [window.tabs.tabText(i) for i in range(window.tabs.count())]
    assert "Planar arrays" in labels
    assert "Linear arrays" in labels, "the linear tab should say which it is"


def test_planar_tab_reports_the_design_sidelobe_level(window):
    tab = window.planar
    tab.lattice.setCurrentText("rectangular")
    tab.taper.setCurrentText("chebyshev")
    for level in (-25.0, -35.0):
        tab.sll.setValue(level)
        rows = _rows(tab.summary_table)
        for key in ("sidelobe, scan plane", "sidelobe, cross plane"):
            assert abs(float(rows[key].split()[0]) - level) < 0.1, (
                f"{key} read {rows[key]} for a {level} dB design")


def test_planar_tab_warns_about_grating_lobes(window):
    tab = window.planar
    tab.lattice.setCurrentText("rectangular")
    tab.scan.setValue(45.0)
    tab.spacing.setValue(0.5)
    assert "grating-lobe limit" not in tab.warning.text()
    tab.spacing.setValue(0.9)
    assert "grating lobe is in real space" in tab.warning.text()
    tab.spacing.setValue(0.5)


def test_planar_tab_says_a_triangular_lattice_cannot_be_tapered(window):
    """Silently ignoring the taper would be worse than saying so."""
    tab = window.planar
    tab.taper.setCurrentText("chebyshev")
    tab.lattice.setCurrentText("triangular")
    assert not tab.taper.isEnabled()
    assert "not separable" in tab.warning.text()
    rows = _rows(tab.summary_table)
    assert float(rows["taper efficiency"]) == pytest.approx(1.0, abs=1e-6)
    tab.lattice.setCurrentText("rectangular")


def test_planar_tab_scan_broadens_only_the_scan_plane(window):
    tab = window.planar
    tab.lattice.setCurrentText("rectangular")
    tab.taper.setCurrentText("uniform")
    tab.spacing.setValue(0.5)
    tab.scan.setValue(0.0)
    base = _rows(tab.summary_table)
    tab.scan.setValue(45.0)
    scanned = _rows(tab.summary_table)
    def deg(rows, key):
        return float(rows[key].rstrip("°"))
    assert deg(scanned, "beamwidth, scan plane") > deg(base, "beamwidth, scan plane") * 1.3
    assert deg(scanned, "beamwidth, cross plane") == pytest.approx(
        deg(base, "beamwidth, cross plane"), rel=2e-3)
    tab.scan.setValue(0.0)


def test_planar_tab_ground_plane_box_doubles_the_directivity(window):
    tab = window.planar
    tab.ground.setChecked(False)
    open_space = float(_rows(tab.summary_table)["directivity"].split()[0])
    assert "mirror beam" in tab.warning.text()
    tab.ground.setChecked(True)
    backed = float(_rows(tab.summary_table)["directivity"].split()[0])
    assert backed - open_space == pytest.approx(3.0103, abs=0.01)


def test_planar_tab_survives_a_one_element_array(window):
    """Degenerate settings must not raise; the tab has to keep working."""
    tab = window.planar
    tab.nx.setValue(1)
    tab.ny.setValue(1)
    rows = _rows(tab.summary_table)
    assert rows["elements"] == "1"
    tab.nx.setValue(12)
    tab.ny.setValue(12)


# ------------------------------------------------------------ the redesign: pictures and the gallery

def test_every_archetype_has_its_own_drawing(registry):
    """No archetype falls back to its family's generic picture."""
    from otahub.gui.drawings import DRAWINGS
    missing = [a.key for a in registry if a.key not in DRAWINGS]
    assert not missing, f"no drawing for {missing}"


def test_every_drawing_paints_from_its_design(qapp, registry):
    """Each figure renders from its own default design, puts copper or metal on
    the canvas, and survives a design that produced nothing."""
    from PySide6.QtCore import QRectF
    from PySide6.QtGui import QColor, QImage, QPainter
    from otahub.gui import drawings
    from otahub.gui.app import default_values, design_values
    blank = []
    for a in registry:
        values = design_values(a.synthesize(**default_values(a)))
        for vals in (values, {}):
            img = QImage(360, 260, QImage.Format.Format_ARGB32)
            img.fill(QColor("#ffffff"))
            p = QPainter(img)
            drawings.paint(p, QRectF(0, 0, 360, 260), a.key, a.family, vals, grid=False)
            p.end()
            inked = sum(1 for x in range(0, 360, 4) for y in range(0, 260, 4) if img.pixelColor(x, y) != QColor("#ffffff"))
            if inked < 20:
                blank.append((a.key, bool(vals)))
    assert not blank, f"drawings that painted nothing: {blank}"


def test_the_gallery_has_one_card_per_archetype(window, registry):
    assert set(window.catalogue.cards) == set(registry.keys)


def test_the_family_chips_filter_the_gallery(window, registry):
    tab = window.catalogue
    tab._on_family("horn")
    shown = {k for k, c in tab.cards.items() if not c.isHidden()}
    assert shown == {a.key for a in registry.by_family("horn")}
    tab._on_family("all")
    assert all(not c.isHidden() for c in tab.cards.values())


def test_search_and_family_combine(window):
    tab = window.catalogue
    tab._on_family("patch")
    tab.search.setText("annular")
    shown = {k for k, c in tab.cards.items() if not c.isHidden()}
    assert shown == {"annular_ring_patch"}
    tab.search.setText("")
    tab._on_family("all")


def test_the_form_understands_units(window):
    tab = window.catalogue
    tab.select_key("rectangular_patch")
    tab._fields["f0"].setText("2.4 GHz")
    tab._fields["h"].setText("1.6 mm")
    values = tab.collect()
    assert values["f0"] == pytest.approx(2.4e9) and values["h"] == pytest.approx(1.6e-3)


def test_quantities_round_trip_through_the_form():
    from otahub.gui.models import format_input, parse_quantity
    for value, unit in ((2.4e9, "Hz"), (5.47723e9, "Hz"), (0.0016, "m"), (8.1e-7, "m"), (4.4, "-"), (50.0, "ohm")):
        assert parse_quantity(format_input(value, unit), unit) == pytest.approx(value, rel=1e-6)
    for bad in ("abc", "2.4 furlong", "1..2"):
        with pytest.raises(ValueError, match="not a number"):
            parse_quantity(bad, "Hz")


def test_headline_figures_never_show_a_loss_resistance_as_the_impedance():
    from otahub.gui.models import key_figures
    figs = dict((label, key) for label, _, key in key_figures(
        {"loss_resistance_ohm": 0.2, "edge_resistance_ohm": 300.0, "inset_resistance_ohm": 50.0,
         "directivity_dbi": 7.0}, {"inset_resistance_ohm": "ohm", "directivity_dbi": "dBi"}))
    assert figs == {"Directivity": "directivity_dbi", "Impedance": "inset_resistance_ohm"}


def test_headline_figures_skip_what_is_outside_the_solved_range():
    from otahub.gui.models import key_figures
    labels = [label for label, _, _ in key_figures({"gain_dbi": float("nan"), "hpbw_e_deg": 30.0},
                                                   {"gain_dbi": "dBi", "hpbw_e_deg": "deg"})]
    assert labels == ["Beamwidth"]


def test_the_drawing_follows_the_design(window):
    tab = window.catalogue
    tab.select_key("half_wave_dipole")
    tab._fields["f0"].setText("300 MHz")
    tab._synthesise()
    long_arm = tab.drawing.values["L"]
    tab._fields["f0"].setText("600 MHz")
    tab._synthesise()
    assert tab.drawing.values["L"] == pytest.approx(long_arm / 2, rel=1e-6)


def test_a_bad_entry_is_marked_on_the_page_not_in_a_dialog(window):
    tab = window.catalogue
    tab.select_key("half_wave_dipole")
    tab._fields["f0"].setText("fast")
    tab._synthesise()
    assert tab._fields["f0"].property("invalid") == "true"
    assert not tab.banner.isHidden() and "f0" in tab.banner_text.text()


def test_the_sweep_offers_readable_quantities(window):
    tab = window.catalogue
    tab.select_key("half_wave_dipole")
    assert tab.metric_picker.currentData() == "directivity_dbi"
    assert tab.metric_picker.currentText() == "Directivity"
    assert tab.sweep_canvas.axes.lines, "the sweep should have drawn a curve"


# ------------------------------------------------------------ found by review, fixed

def _labels(key, registry, values=None):
    from PySide6.QtCore import QRectF
    from PySide6.QtGui import QColor, QImage, QPainter
    from otahub.gui import drawings
    from otahub.gui.app import design_values
    a = registry[key]
    if values is None:
        values = design_values(a.synthesize(**a.spec.known_cases[0].given))
    img = QImage(520, 400, QImage.Format.Format_ARGB32)
    img.fill(QColor("#ffffff"))
    p = QPainter(img)
    try:
        return drawings.paint(p, QRectF(0, 0, 520, 400), key, a.family, values)
    finally:
        p.end()


@pytest.mark.parametrize("key,label", [
    ("lpda", "28 elements"),                                  # was 16, at a scale factor of its own
    ("normal_mode_helix", "40 mm, 20 turns"),                 # was four turns
    ("axial_mode_helix", "708 mm, 10 turns"),
    ("waveguide_slot_array_resonant", "12 slots"),            # was always six
    ("waveguide_slot_array_travelling_wave", "20 slots"),
    ("cylindrical_parabolic", "W = 1 m"),                     # was D = 1 m from a default
    ("cassegrain", "F = 1.05 m"),                             # f/D times D, not a fixed proportion
    ("conical_horn", "slant 300 mm"),                         # along the wall from the apex
    ("truncated_corner_cp_patch", "cut 3.97 mm"),             # a leg, not the diagonal
    ("vivaldi_tsa", "aperture 150 mm"),
])
def test_each_drawing_labels_its_own_designs_numbers(qapp, registry, key, label):
    """The drawings replaced counts, angles and sizes with fixed ones; the labels
    they paint must be the design's."""
    assert label in _labels(key, registry)


def _dim(record, prefix):
    x1, y1, x2, y2, _ = next(d for d in record.dims if d[4].startswith(prefix))
    return x1, y1, x2, y2, math.hypot(x2 - x1, y2 - y1)


def test_the_drawings_geometry_is_the_designs_not_just_its_labels(qapp, registry):
    """The labels were already right where the lines were wrong (a review swapped
    the old renderers back in and the label tests passed), so check where the
    lines go: the CP cut spans one leg, the conical slant runs apex to rim at its
    true length, the Vivaldi's aperture line spans W_ap, the LPDA draws every
    element."""
    from otahub.gui.app import design_values
    cp = registry["truncated_corner_cp_patch"]
    d = cp.synthesize(**cp.spec.known_cases[0].given)
    x1, y1, x2, y2, n = _dim(_labels("truncated_corner_cp_patch", registry), "cut ")
    assert y1 == pytest.approx(y2) and n == pytest.approx(d.get("c_trunc"), rel=1e-9)
    horn = registry["conical_horn"]
    d = horn.synthesize(**horn.spec.known_cases[0].given)
    x1, y1, x2, y2, n = _dim(_labels("conical_horn", registry), "slant ")
    assert n == pytest.approx(d.get("L"), rel=1e-9) and y1 == pytest.approx(0.0)
    long = registry["vivaldi_tsa"].synthesize(f_low=1e9, Lax=2.0)     # a board wider than its aperture
    rec = _labels("vivaldi_tsa", registry, design_values(long))
    assert _dim(rec, "aperture ")[4] == pytest.approx(long.get("W_ap"), rel=1e-9)
    big = registry["lpda"].synthesize(f_low=10e6, f_high=10e9, tau=0.92)
    rec = _labels("lpda", registry, design_values(big))
    assert len(rec.wires) == 2 * int(big.get("N_elements")) == 178


def test_a_three_level_zone_plate_steps_every_two_thirds_of_a_zone(qapp, registry):
    """levels // 2 lost the fraction: four bands from 39.04 mm where the model has
    six from 31.79 mm, r_m = sqrt(m lambda F + (m lambda / 2)^2) at m = 2j/3."""
    from otahub.gui.app import design_values
    d = registry["fresnel_zone_plate"].synthesize(f0=30e9, F=0.15, M=4, phase_levels=3)
    rec = _labels("fresnel_zone_plate", registry, design_values(d))
    lam = 2.99792458e8 / 30e9
    want = [math.sqrt(m * lam * 0.15 + (m * lam / 2) ** 2) for m in (2 * j / 3 for j in range(1, 7))]
    assert sorted(rec.circles) == pytest.approx(want, rel=1e-12)
    assert min(rec.circles) == pytest.approx(0.0317869, rel=1e-5)


@pytest.mark.parametrize("key,kind", [("cassegrain", "cass"), ("gregorian_dual_reflector", "greg")])
def test_the_subreflector_is_the_conic_its_magnification_makes(qapp, registry, key, kind):
    """It was a fixed parabola at a fixed place whatever the magnification. Now a
    hyperbola (Cassegrain) or ellipse (Gregorian) on the two foci, whose
    eccentricity is the magnification's: (M+1)/(M-1), or (M-1)/(M+1)."""
    from otahub.gui.app import design_values
    a = registry[key]
    shapes = []
    for M in (4.0, 8.0):
        g = dict(a.spec.known_cases[0].given, magnification=M)
        dual = _labels(key, registry, design_values(a.synthesize(**g))).meta["dual"]
        e = (M + 1) / (M - 1) if kind == "cass" else (M - 1) / (M + 1)
        assert dual["e"] == pytest.approx(e, rel=1e-12)
        f1, f2 = dual["f1"], dual["f2"]
        focal = [math.hypot(x - f2, y) + (-1 if kind == "cass" else 1) * math.hypot(x - f1, y)
                 for x, y in dual["sub"]]
        assert max(focal) - min(focal) < 1e-12
        assert max(abs(y) for _, y in dual["sub"]) == pytest.approx(g["Ds"] / 2 if "Ds" in g else dual["rs"])
        shapes.append(dual["sub"])
    assert shapes[0] != shapes[1]


def test_a_drawing_never_labels_a_number_the_design_did_not_compute(qapp, registry):
    """A dipole with no length was labelled "L = 1 m". With no design at all the only
    numbers left are the archetypes' own fixed angles."""
    fixed = {"90°", "60°", "90° phasing line"}
    for a in registry:
        numbers = [s for s in _labels(a.key, registry, {}) if re.search(r"\d", s) and s not in fixed]
        assert not numbers, (a.key, numbers)


def test_every_input_a_design_asks_for_is_on_the_form(window, registry):
    """A small square loop needs its turn count and wire radius, and the form had no
    field for either: geometry the synthesis reads but never produces is an input."""
    tab = window.catalogue
    for a in registry:
        tab.select_key(a.key)
        missing = set(tab._design.missing_requirements()) if tab._design else set()
        assert missing <= set(tab._fields), (a.key, missing - set(tab._fields))
    tab.select_key("small_square_loop")
    assert {"N", "b"} <= set(tab._fields)


def test_every_page_opens_on_a_complete_design(window, registry):
    """Seven loop pages opened with a blank wire radius (the ferrite with no rod
    diameter or coil length) and so no impedance at all. They open complete now:
    from the known design that fills most of the page, and otherwise from the
    middle of the spec's typical range, worked out on the design."""
    tab = window.catalogue
    incomplete = {}
    for a in registry:
        tab.select_key(a.key)
        missing = tab._design.missing_requirements() if tab._design else ["no design"]
        if missing:
            incomplete[a.key] = missing
    assert not incomplete, incomplete
    tab.select_key("small_square_loop")
    s = tab._design.get("s")
    assert 0.001 * s < tab._design.get("b") < 0.05 * s


def test_loading_a_known_design_resets_what_it_leaves_out(window):
    """A setting changed by hand stayed when a published design was loaded, so the
    page showed a design that was not the published one."""
    tab = window.catalogue
    tab.select_key("cylindrical_parabolic")
    tab._fields["f_over_W"].setText("0.8")
    tab._load_example(1)
    assert tab._design.get("F") == pytest.approx(0.4)


def test_leaving_the_design_page_stops_a_pending_recalculation(window):
    tab = window.catalogue
    tab.select_key("half_wave_dipole")
    tab._timer.start()
    tab.show_gallery()
    assert not tab._timer.isActive()
    tab.select_key("half_wave_dipole")
    tab._timer.start()
    window._go(2)
    assert not tab._timer.isActive()


def test_the_monopole_pattern_radiates_only_above_its_ground(registry):
    """Drawn as the whole image dipole it read 2.15 dBi with a lower half; the spec
    says 5.16 dBi."""
    from otahub.gui.app import PATTERN_SOURCES
    pattern, _ = PATTERN_SOURCES["quarter_wave_monopole"]()
    below = pattern.U[pattern.theta > math.pi / 2]
    assert not below.any()
    spec = registry["quarter_wave_monopole"].synthesize(f0=100e6).metrics
    d = next(v for k, v in spec.items() if k.startswith("directivity") and k.endswith("dbi"))
    assert pattern.directivity_dbi() == pytest.approx(d, abs=0.02)


def test_the_planar_warning_tests_the_steering_direction(window):
    """0.6 lambda at 45 deg in the phi = 45 plane has no grating lobe, though it is
    past the all-azimuth limit; the page said a lobe was in real space."""
    tab = window.planar
    tab.lattice.setCurrentText("rectangular")
    tab.spacing.setValue(0.6)
    tab.scan.setValue(45)
    tab.scan_phi.setValue(45)
    tab.refresh()
    assert "A grating lobe is in real space" not in tab.warning.text()
    assert "in every plane" in tab.warning.text()
    tab.scan_phi.setValue(0)
    tab.refresh()
    assert "A grating lobe is in real space" in tab.warning.text()

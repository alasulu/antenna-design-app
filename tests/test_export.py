"""Exporter tests.

The point of these is not that the scripts look plausible - it is that they
are structurally complete, that units are right, and that an archetype we
cannot build says so instead of emitting a half-model.
"""
import re

import math

import pytest

from otahub.export import build
from otahub.export import cst, hfss
from otahub.export.base import BUILDERS, Brick, Cylinder, DiscretePort

CASES = {
    "half_wave_dipole": {"f0": 300e6, "aw": 1e-3},
    "short_dipole": {"f0": 300e6, "L_over_lambda": 0.01},
    "resonant_dipole": {"f0": 300e6, "aw": 1e-3},
    "quarter_wave_monopole": {"f0": 1e9, "aw": 1e-3},
    "rectangular_patch": {"f0": 2.4e9, "eps_r": 4.4, "h": 1.6e-3},
    "rectangular_patch_inset": {"f0": 2.4e9, "eps_r": 4.4, "h": 1.6e-3, "Z_target": 50.0},
    "circular_patch": {"f0": 10e9, "eps_r": 2.2, "h": 1.588e-3},
    "open_ended_waveguide": {"f0": 10e9, "a_wg": 0.02286, "b_wg": 0.01016},
    "folded_dipole": {"f0": 300e6, "N": 2, "aw": 1e-3},
    "dipole_over_ground": {"f0": 300e6, "h_over_lambda": 0.25, "aw": 1e-3},
    "turnstile_dipole": {"f0": 300e6, "aw": 1e-3},
    "quarter_wave_shorted_patch": {"f0": 2.4e9, "eps_r": 4.4, "h": 1.6e-3},
    "rectangular_dra": {"f0": 10e9, "eps_r": 10.0},
    "cylindrical_dra": {"f0": 10e9, "eps_r": 10.0},
    "hemispherical_dra": {"f0": 10e9, "eps_r": 10.0},
    "conical_monopole": {"f_low": 1e9, "cone_half_angle_deg": 47.0},
    "biconical": {"f0": 1e9},
    "discone": {"f_low": 100e6},
    "small_circular_loop": {"f0": 10e6, "C_over_lambda": 0.1},
    "one_wavelength_circular_loop": {"f0": 300e6},
    "small_square_loop": {"f0": 10e6, "P_over_lambda": 0.1},
    "quad_loop_square": {"f0": 144e6},
    "alford_loop": {"f0": 300e6, "P_over_lambda": 0.5},
    "halo_loop": {"f0": 144e6},
    "half_wave_slot": {"f0": 300e6},
    "folded_slot": {"f0": 300e6, "N": 2},
    "cavity_backed_slot": {"f0": 2.4e9},
    "waveguide_longitudinal_slot": {"f0": 10e9, "a_wg": 0.02286,
                                    "b_wg": 0.01016, "x1": 0.003},
    "waveguide_slot_array_resonant": {"f0": 10e9, "a_wg": 0.02286,
                                      "b_wg": 0.01016, "N": 12},
    "waveguide_slot_array_travelling_wave": {"f0": 10e9, "a_wg": 0.02286,
                                             "b_wg": 0.01016},
    "long_wire_travelling": {"f0": 300e6, "L_over_lambda": 4.0},
    "leaky_wave_line_source": {"f0": 10e9},
    "planar_monopole_rectangular": {"f_low": 1.5e9},
    "planar_monopole_circular": {"f_low": 1.5e9},
    "annular_ring_patch": {"f0": 2e9, "eps_r": 2.2, "h": 1.6e-3},
    "pifa": {"f0": 2.4e9, "h": 0.006},
    "stacked_patch": {"f0": 2.4e9, "eps_r": 2.2, "h": 1.6e-3},
    "fresnel_zone_plate": {"f0": 30e9, "F": 0.15, "M": 4},
    "luneburg_lens": {"f0": 30e9, "D": 0.3},
    "dipole_arbitrary_length": {"f0": 300e6, "L_over_lambda": 0.75, "aw": 1e-4},
    "top_loaded_monopole": {"f0": 10e6, "h_over_lambda": 0.05, "beta_top": 0.6},
    "inductively_loaded_monopole": {"f0": 10e6, "h_over_lambda": 0.05},
    "multiturn_small_loop": {"f0": 10e6, "C_over_lambda": 0.1, "Rin_target": 50.0},
}


@pytest.fixture(params=sorted(CASES))
def model(request, registry):
    key = request.param
    return build(registry[key].synthesize(**CASES[key]))


# -------------------------------------------------------------- geometry IR

def test_every_registered_builder_is_exercised_by_a_test_case():
    assert set(BUILDERS) <= set(CASES), f"untested builders: {set(BUILDERS) - set(CASES)}"


def test_builders_produce_solids(model):
    assert model.built_geometry
    assert model.solids


def test_builders_carry_units_for_every_parameter(model):
    missing = [n for n in model.parameters if n not in model.units]
    assert not missing, f"no declared unit for {missing}"


def test_dipole_arms_are_symmetric_about_the_feed(registry):
    m = build(registry["half_wave_dipole"].synthesize(f0=300e6, aw=1e-3))
    upper = next(s for s in m.solids if s.name == "arm_upper")
    lower = next(s for s in m.solids if s.name == "arm_lower")
    assert upper.span[1] == pytest.approx(-lower.span[0])
    assert upper.span[0] == pytest.approx(-lower.span[1])


def test_dipole_total_length_matches_the_synthesised_value(registry):
    design = registry["half_wave_dipole"].synthesize(f0=300e6, aw=1e-3)
    m = build(design)
    upper = next(s for s in m.solids if s.name == "arm_upper")
    lower = next(s for s in m.solids if s.name == "arm_lower")
    tip_to_tip = upper.span[1] - lower.span[0]
    assert tip_to_tip == pytest.approx(design.get("L"), rel=1e-9)


def test_patch_ground_sits_under_the_substrate(registry):
    m = build(registry["rectangular_patch"].synthesize(f0=2.4e9, eps_r=4.4, h=1.6e-3))
    ground = next(s for s in m.solids if s.name == "ground")
    patch = next(s for s in m.solids if s.name == "patch")
    assert ground.z == (0.0, 0.0)
    assert patch.z[0] == pytest.approx(1.6e-3)


def test_unbuildable_archetype_reports_instead_of_faking(registry):
    """The honesty guarantee: no half-model that looks complete."""
    m = build(registry["pyramidal_horn"].synthesize(f0=10e9, G_target=20.0))
    assert not m.built_geometry
    assert not m.solids
    assert m.parameters                       # parameters still exported
    assert any("NO SOLID GEOMETRY" in n for n in m.notes)


def test_every_archetype_exports_without_raising(registry):
    """Whatever we cannot build must still degrade to a parameters-only model."""
    import math
    for archetype in registry:
        lo, hi = archetype.spec.freq_range_hz
        f0 = math.sqrt(max(lo, 1e5) * min(hi if math.isfinite(hi) else 1e11, 1e11))
        reqs = {"f0": f0}
        for p in archetype.spec.parameters:
            if p.role in ("requirement", "assumption", "material") and p.symbol != "f0":
                try:
                    reqs[p.symbol] = float(p.typical)
                except (TypeError, ValueError):
                    continue
        m = build(archetype.synthesize(**reqs))
        assert cst.render(m)
        assert hfss.render(m)


# ------------------------------------------------------------------ backends

def test_cst_macro_is_structurally_complete(model):
    text = cst.render(model)
    assert text.count("Sub Main") == 1
    assert text.count("End Sub") == 1
    assert text.index("Sub Main") < text.index("End Sub")
    assert text.count("With ") == text.count("End With")


def test_cst_creates_every_solid(model):
    text = cst.render(model)
    for solid in model.solids:
        assert f'.Name "{solid.name}"' in text


def test_cst_declares_a_discrete_port_when_one_exists(model):
    text = cst.render(model)
    if model.ports:
        assert "With DiscretePort" in text
        assert ".Create" in text


def test_hfss_script_is_structurally_complete(model):
    text = hfss.render(model)
    assert "ScriptEnv.Initialize" in text
    assert "oEditor = oDesign.SetActiveEditor" in text
    assert text.count("(") == text.count(")")
    assert text.count("[") == text.count("]")


def test_hfss_script_is_valid_python_syntax(model):
    """It is IronPython, but the syntax must still parse."""
    import ast
    body = hfss.render(model)
    body = body.replace("import ScriptEnv", "ScriptEnv = None")
    ast.parse(body)


def test_hfss_creates_every_solid(model):
    text = hfss.render(model)
    for solid in model.solids:
        assert f'"Name:=", "{solid.name}"' in text


# ---------------------------------------------------------------------- units

NON_LENGTH_UNITS = {"ohm", "S", "-", "", "Hz", "F", "H", "rad/m", "dBi", "dB"}


def test_no_non_length_parameter_is_written_as_a_length(model):
    """Guards the bug where a 319 ohm resistance exported as '319105 mm'."""
    text = cst.render(model)
    for name, unit in model.units.items():
        if name not in model.parameters or unit == "m":
            continue
        assert f'"{name}_mm"' not in text, f"{name} [{unit}] exported as a length"


def test_lengths_are_converted_to_millimetres(registry):
    design = registry["half_wave_dipole"].synthesize(f0=300e6, aw=1e-3)
    m = build(design)
    text = cst.render(m)
    match = re.search(r'StoreParameter "L_mm", ([0-9.]+)', text)
    assert match
    assert float(match.group(1)) == pytest.approx(design.get("L") * 1e3, rel=1e-6)


def test_frequency_is_exported_in_ghz(registry):
    m = build(registry["half_wave_dipole"].synthesize(f0=2.4e9, aw=1e-3))
    text = cst.render(m)
    assert re.search(r'StoreParameter "f0_GHz", 2\.4', text)
    for backend_text in (text, hfss.render(m)):
        assert "2400000000" not in backend_text


def test_impedance_keeps_its_own_magnitude(registry):
    design = registry["rectangular_patch_inset"].synthesize(
        f0=2.4e9, eps_r=4.4, h=1.6e-3, Z_target=50.0)
    text = cst.render(build(design))
    match = re.search(r'StoreParameter "Rin0", ([0-9.]+)', text)
    assert match
    assert float(match.group(1)) == pytest.approx(design.get("Rin0"), rel=1e-6)


def test_notes_are_carried_into_both_backends(model):
    assert all(any(w in cst.render(model) for w in note.split()[:3])
               for note in model.notes[:1]) or not model.notes
    if model.notes:
        assert "NOTE:" in cst.render(model)
        assert "NOTE:" in hfss.render(model)


# ------------------------------------------------- cones, spheres, booleans

def test_boolean_operations_only_reference_solids_that_exist(model):
    """A subtract naming a solid that was never created fails silently in both
    tools: CST reports a missing object and HFSS raises inside the script."""
    from otahub.export.base import Subtract

    names = {s.name for s in model.solids}
    for op in model.operations:
        if isinstance(op, Subtract):
            assert op.target in names, f"subtract target {op.target!r} is not a solid"
            for tool in op.tools:
                assert tool in names, f"subtract tool {tool!r} is not a solid"


def test_boolean_operations_are_emitted_after_every_solid(model):
    """Order matters: both backends operate on named objects, so a subtract
    emitted before its operands is a runtime error in the simulator."""
    from otahub.export.base import Subtract

    subs = [op for op in model.operations if isinstance(op, Subtract)]
    if not subs:
        pytest.skip("no boolean operations in this model")
    # locate each solid's CREATION, not its name's first mention: the header and the
    # parameter list mention names too, and let a solid created after its boolean pass
    for text, marker, created in ((cst.render(model), "Solid.Subtract", '        .Name "{}"'),
                                  (hfss.render(model), "oEditor.Subtract", '"Name:=", "{}"')):
        first_op = text.index(marker)
        for solid in model.solids:
            assert text.index(created.format(solid.name)) < first_op, (
                f"{solid.name} is created after the boolean that uses it")


def test_cone_radii_reach_both_backends(registry):
    from otahub.export.base import Cone

    model = build(registry["discone"].synthesize(f_low=100e6))
    cones = [s for s in model.solids if isinstance(s, Cone)]
    assert cones, "the discone must contain a cone"
    cone = cones[0]
    # truncated at the feed: the spec's band edges were solved with that top face
    assert cone.radius_start > cone.radius_end > 0
    vba, py = cst.render(model), hfss.render(model)
    assert "With Cone" in vba and "Bottomradius" in vba and "Topradius" in vba
    assert "CreateCone" in py and "BottomRadius:=" in py and "TopRadius:=" in py


def test_sphere_radius_reaches_both_backends(registry):
    from otahub.export.base import Sphere

    model = build(registry["hemispherical_dra"].synthesize(f0=10e9, eps_r=10.0))
    spheres = [s for s in model.solids if isinstance(s, Sphere)]
    assert spheres, "the hemispherical DRA must contain a sphere"
    radius_mm = spheres[0].radius * 1e3
    assert f"{radius_mm:.6f}" in cst.render(model)
    assert "CreateSphere" in hfss.render(model)


def test_biconical_cone_has_the_specs_half_angle(registry):
    """theta_h is the half angle from the axis. The builder once halved it
    again and drew a 15 degree cone for a 30 degree design."""
    from otahub.export.base import Cone

    design = registry["biconical"].synthesize(f0=1e9, theta_h=math.radians(30.0))
    cone = [s for s in build(design).solids if isinstance(s, Cone) and s.name == "cone_upper"][0]
    height = cone.span[1] - cone.span[0]
    flare = math.atan((cone.radius_end - cone.radius_start) / height)
    assert math.degrees(flare) == pytest.approx(30.0, abs=1e-6)


def test_biconical_cones_are_mirror_images(registry):
    from otahub.export.base import Cone

    model = build(registry["biconical"].synthesize(f0=1e9))
    upper, lower = [s for s in model.solids if isinstance(s, Cone)]
    assert upper.radius_end == pytest.approx(lower.radius_start)
    assert upper.span[1] - upper.span[0] == pytest.approx(lower.span[1] - lower.span[0])
    assert upper.span[0] == pytest.approx(-lower.span[1])


def test_turnstile_exports_two_ports(registry):
    model = build(registry["turnstile_dipole"].synthesize(f0=300e6, aw=1e-3))
    assert len(model.ports) == 2, "a turnstile needs both dipoles driven"
    py = hfss.render(model)
    assert py.count("AssignLumpedPort") == 2


def test_dra_resonator_sits_on_the_ground_plane(registry):
    """The DRA family all assume a ground plane at z = 0; a resonator floating
    above it or buried below models a different antenna."""
    for key, given in (("rectangular_dra", {"f0": 10e9, "eps_r": 10.0}),
                       ("cylindrical_dra", {"f0": 10e9, "eps_r": 10.0})):
        model = build(registry[key].synthesize(**given))
        res = [s for s in model.solids if s.name == "resonator"][0]
        span = res.z if hasattr(res, "z") else res.span
        assert span[0] == pytest.approx(0.0), f"{key} resonator starts at {span[0]}"
        assert span[1] > 0


# ------------------------------------------------------------------- loops

def test_torus_reaches_both_backends_with_the_right_radii(registry):
    """CST wants inner and outer radii measured from the axis; HFSS wants major
    and minor. Feeding either one the other's numbers makes a torus of the
    wrong size, and nothing in the file looks wrong."""
    from otahub.export.base import Torus

    model = build(registry["one_wavelength_circular_loop"].synthesize(f0=300e6))
    tori = [s for s in model.solids if isinstance(s, Torus)]
    assert tori, "a circular loop must contain a torus"
    t = tori[0]
    vba, py = cst.render(model), hfss.render(model)
    assert "With Torus" in vba and "CreateTorus" in py
    # CST: outer and inner, derived
    assert f"{(t.major_radius + t.minor_radius)*1e3:.6f}" in vba
    assert f"{(t.major_radius - t.minor_radius)*1e3:.6f}" in vba
    # HFSS: major and minor, as given
    assert f'"MajorRadius:=", "{t.major_radius*1e3:.6f}' in py
    assert f'"MinorRadius:=", "{t.minor_radius*1e3:.6f}' in py


def test_loop_feed_gap_is_cut_not_merely_drawn(registry):
    """A port across an unbroken ring shorts itself out. The gap has to be a
    real boolean subtraction."""
    from otahub.export.base import Subtract

    for key, given in (("small_circular_loop", {"f0": 10e6, "C_over_lambda": 0.1}),
                       ("halo_loop", {"f0": 144e6})):
        model = build(registry[key].synthesize(**given))
        subs = [op for op in model.operations if isinstance(op, Subtract)]
        assert subs, f"{key} must cut its feed gap"
        names = {s.name for s in model.solids}
        assert subs[0].target in names
        assert all(t in names for t in subs[0].tools)


def test_square_loop_sides_close_the_perimeter(registry):
    """Four sides, one of them split for the feed - five cylinders whose spans
    add up to the perimeter less the gap."""
    from otahub.export.base import Cylinder

    design = registry["quad_loop_square"].synthesize(f0=144e6)
    model = build(design)
    cyls = [s for s in model.solids if isinstance(s, Cylinder)]
    assert len(cyls) == 5, "three whole sides plus a split fourth"
    total = sum(abs(c.span[1] - c.span[0]) for c in cyls)
    side = design.get("s")
    gap = 4 * side - total
    assert 0 < gap < side / 10, f"gap {gap} is not a small fraction of a side"


def test_halo_gap_comes_from_the_spec_not_from_the_builder(registry):
    """Every other loop's feed gap is invented by the exporter. The halo's TIP gap
    is a design parameter, because its capacitance sets the resonance - and it is
    left open: the halo is fed at the middle of the conductor, opposite the tips,
    as the verified MoM model drives it. The port used to sit across the tip gap."""
    design = registry["halo_loop"].synthesize(f0=144e6)
    model = build(design)
    cut = next(s for s in model.solids if s.name == "tip_gap_cut")
    assert cut.y[1] - cut.y[0] == pytest.approx(design.get("g"), rel=1e-9)
    port = model.ports[0]
    radius = design.get("Dm") / 2
    assert port.start[0] == pytest.approx(-radius) and port.end[0] == pytest.approx(-radius)
    assert min(s.x[1] for s in model.solids if s.name == "tip_gap_cut") > 0   # the tips at +x


def test_loop_builders_say_what_they_had_to_invent(registry):
    """Feed gaps, wire radii and missing corner capacitors all change the
    answer; a model that stays silent about them is misleading."""
    for key, given in (("small_circular_loop", {"f0": 10e6, "C_over_lambda": 0.1}),
                       ("small_square_loop", {"f0": 10e6, "P_over_lambda": 0.1}),
                       ("alford_loop", {"f0": 300e6, "P_over_lambda": 0.5})):
        model = build(registry[key].synthesize(**given))
        joined = " ".join(model.notes).lower()
        assert "gap" in joined, f"{key} does not mention its feed gap"


# -------------------------------------------------------------------- slots

def test_a_slot_is_cut_from_its_ground_plane(registry):
    """A slot drawn as a separate solid is not a slot, it is a plate with a bar
    lying on it. Every slot archetype must subtract."""
    from otahub.export.base import Subtract

    for key, given in (("half_wave_slot", {"f0": 300e6}),
                       ("folded_slot", {"f0": 300e6, "N": 2}),
                       ("cavity_backed_slot", {"f0": 2.4e9})):
        model = build(registry[key].synthesize(**given))
        subs = [op for op in model.operations if isinstance(op, Subtract)]
        assert subs, f"{key} never cuts its slot"
        assert subs[0].target == "ground_plane"


def test_folded_slot_cuts_one_slot_per_conductor(registry):
    for n in (1, 2, 3):
        model = build(registry["folded_slot"].synthesize(f0=300e6, N=n))
        cuts = [s for s in model.solids if s.name.startswith("slot_cut")]
        assert len(cuts) == n
        assert len(model.operations[0].tools) == n


def test_slot_array_alternates_its_offsets(registry):
    """The alternation is the whole mechanism: it undoes the 180 degrees of
    propagation phase between slots half a guide wavelength apart. Building
    them all on one side gives two beams off broadside instead of one on it."""
    design = registry["waveguide_slot_array_resonant"].synthesize(
        f0=10e9, a_wg=0.02286, b_wg=0.01016, N=12)
    model = build(design)
    cuts = sorted((s for s in model.solids if s.name.startswith("slot_cut")),
                  key=lambda s: s.z[0])
    assert len(cuts) == 12
    centres = [(c.x[0] + c.x[1]) / 2 for c in cuts]
    signs = [1 if c > 0 else -1 for c in centres]
    assert all(a != b for a, b in zip(signs, signs[1:])), "offsets must alternate"
    assert all(abs(abs(c) - design.get("offset")) < 1e-12 for c in centres)


def test_slot_array_spacing_matches_the_spec(registry):
    for key, given, want in (
            ("waveguide_slot_array_resonant",
             {"f0": 10e9, "a_wg": 0.02286, "b_wg": 0.01016, "N": 12}, "spacing"),
            ("waveguide_slot_array_travelling_wave",
             {"f0": 10e9, "a_wg": 0.02286, "b_wg": 0.01016}, "spacing")):
        design = registry[key].synthesize(**given)
        model = build(design)
        cuts = sorted((s for s in model.solids if s.name.startswith("slot_cut")),
                      key=lambda s: s.z[0])
        gaps = [b.z[0] - a.z[0] for a, b in zip(cuts, cuts[1:])]
        assert all(g == pytest.approx(design.get(want), rel=1e-9) for g in gaps)


def test_slots_are_cut_through_the_broad_wall_not_the_interior(registry):
    """A slot that stops short of the wall's outer face is a blind hole."""
    design = registry["waveguide_longitudinal_slot"].synthesize(
        f0=10e9, a_wg=0.02286, b_wg=0.01016, x1=0.003)
    model = build(design)
    wall = next(s for s in model.solids if s.name == "broad_wall")
    cut = next(s for s in model.solids if s.name.startswith("slot_cut"))
    assert cut.y[0] <= wall.y[0] and cut.y[1] >= wall.y[1], (
        "the cut must span the full wall thickness")


def test_cavity_sits_behind_the_slot_not_in_front(registry):
    design = registry["cavity_backed_slot"].synthesize(f0=2.4e9)
    model = build(design)
    cavity = next(s for s in model.solids if s.name == "cavity")
    assert cavity.z[1] <= 0.0, "the cavity belongs behind the ground plane"
    assert cavity.z[1] - cavity.z[0] == pytest.approx(
        design.get("cavity_depth"), rel=1e-6)


# ------------------------------------------------ terminated and planar shapes

def test_long_wire_exports_a_termination_port_not_just_a_feed(registry):
    """A travelling-wave wire with nothing at the far end is a standing-wave
    wire, and its pattern splits. The termination has to be modelled."""
    model = build(registry["long_wire_travelling"].synthesize(
        f0=300e6, L_over_lambda=4.0))
    assert len(model.ports) == 2, "feed and termination"
    feed, term = model.ports
    assert term.impedance != feed.impedance, (
        "the termination is not a 50 ohm measurement port")
    assert any("termination" in n.lower() for n in model.notes)


def test_leaky_wave_slit_stops_short_of_both_ends(registry):
    """A slit running the full length would cut through the port faces."""
    design = registry["leaky_wave_line_source"].synthesize(f0=10e9)
    model = build(design)
    guide = next(s for s in model.solids if s.name == "guide_interior")
    slit = next(s for s in model.solids if s.name == "slit_cut")
    assert slit.z[0] > guide.z[0] and slit.z[1] < guide.z[1]


@pytest.mark.parametrize("key,given", [
    ("planar_monopole_rectangular", {"f_low": 1.5e9}),
    ("planar_monopole_circular", {"f_low": 1.5e9}),
])
def test_planar_monopoles_stand_above_their_ground_plane(key, given, registry):
    """The feed gap is a design parameter for these - the low-frequency cut-off
    depends on it directly - so the radiator must not touch the plane."""
    design = registry[key].synthesize(**given)
    model = build(design)
    radiator = next(s for s in model.solids
                    if s.name in ("plate", "disc"))
    lowest = radiator.z[0] if hasattr(radiator, "z") else \
        radiator.centre[1] - radiator.radius
    assert lowest == pytest.approx(design.get("p_gap"), rel=1e-6, abs=1e-9)


def test_annular_ring_is_cut_not_drawn_as_two_discs(registry):
    from otahub.export.base import Subtract

    design = registry["annular_ring_patch"].synthesize(
        f0=2e9, eps_r=2.2, h=1.6e-3)
    model = build(design)
    subs = [op for op in model.operations if isinstance(op, Subtract)]
    assert subs and subs[0].target == "ring"
    ring = next(s for s in model.solids if s.name == "ring")
    hole = next(s for s in model.solids if s.name == "ring_hole")
    assert ring.radius == pytest.approx(design.get("b_out"), rel=1e-9)
    assert hole.radius == pytest.approx(design.get("a_in"), rel=1e-9)
    assert hole.radius < ring.radius


def test_stacked_patch_parasitic_sits_above_the_driven_one(registry):
    design = registry["stacked_patch"].synthesize(f0=2.4e9, eps_r=2.2, h=1.6e-3)
    model = build(design)
    driven = next(s for s in model.solids if s.name == "driven_patch")
    para = next(s for s in model.solids if s.name == "parasitic_patch")
    assert para.z[0] > driven.z[0]
    assert para.z[0] - driven.z[0] == pytest.approx(design.get("h2"), rel=1e-9)


def test_pifa_shorting_wall_spans_the_declared_width(registry):
    design = registry["pifa"].synthesize(f0=2.4e9, h=0.006)
    model = build(design)
    wall = next(s for s in model.solids if s.name == "shorting_wall")
    assert wall.x[1] - wall.x[0] == pytest.approx(design.get("Ws"), rel=1e-9)
    assert wall.z[1] == pytest.approx(design.get("h"), rel=1e-9)


def test_low_confidence_archetypes_say_so_in_their_exported_model(registry):
    """A model whose spec is flagged low confidence must carry that warning into
    the file, or it arrives in the solver looking as solid as any other."""
    for key, given in (("pifa", {"f0": 2.4e9, "h": 0.006}),
                       ("stacked_patch", {"f0": 2.4e9, "eps_r": 2.2, "h": 1.6e-3}),
                       ("planar_monopole_circular", {"f_low": 1.5e9}),
                       ("halo_loop", {"f0": 144e6})):
        model = build(registry[key].synthesize(**given))
        joined = " ".join(model.notes).lower()
        assert "confidence" in joined, f"{key} does not carry its low-confidence flag"


# ------------------------------------------------------------ found by review, fixed

def test_every_hfss_port_sheet_holds_its_integration_line(registry):
    """The sheet was normal to Y whatever the feed: an x-directed feed got a sheet
    of zero height, a y-directed one a sheet its line did not lie in. Now a square
    as wide as the gap, in a plane containing the line."""
    from otahub.export.hfss import _port_sheet
    for key in BUILDERS:
        a = registry[key]
        for port in build(a.synthesize(**a.spec.known_cases[0].given)).ports:
            corner, side, normal = _port_sheet(port)
            n = "XYZ".index(normal)
            assert side > 0
            for end in (port.start, port.end):
                assert end[n] == pytest.approx(corner[n], abs=1e-12)          # on the sheet's plane
                for i in {0, 1, 2} - {n}:                                   # and inside it
                    assert corner[i] - 1e-12 <= end[i] <= corner[i] + side + 1e-12, (key, port.name)


def test_hfss_uses_the_documented_sweep_and_port_properties(registry):
    """AEDT's InsertFrequencySweep takes RangeType/RangeStart/RangeEnd/RangeCount
    (StartValue/StopValue/Count define nothing), and a lumped port's impedance is
    the top-level Impedance - RenormImp only renormalises the reported S."""
    text = hfss.render(build(registry["long_wire_travelling"].synthesize(
        **registry["long_wire_travelling"].spec.known_cases[0].given)))
    assert '"RangeType:=", "LinearCount"' in text and '"RangeStart:="' in text and '"RangeCount:="' in text
    assert "StartValue" not in text and '"Count:="' not in text
    assert '"Impedance:=", "600ohm"' in text and '"Impedance:=", "50ohm"' in text


def test_cst_ports_are_numbered_in_turn(registry):
    text = cst.render(build(registry["turnstile_dipole"].synthesize(f0=300e6, aw=0.001)))
    assert re.findall(r'\.PortNumber "(\d+)"', text) == ["1", "2"]


def test_the_loss_tangent_is_the_designs(registry):
    """0.02 was imposed on every dielectric; the specs assume lossless ones, and a
    patch with tan_d carries its own."""
    lossless = build(registry["rectangular_dra"].synthesize(
        **registry["rectangular_dra"].spec.known_cases[0].given))
    assert '.TanD "0"' in cst.render(lossless) and '"dielectric_loss_tangent:=", "0"' in hfss.render(lossless)
    lossy = build(registry["rectangular_patch"].synthesize(f0=2.4e9, eps_r=4.4, h=0.0016, tan_d=0.02))
    assert '.TanD "0.02"' in cst.render(lossy)


def test_a_wideband_design_is_simulated_over_its_own_band(registry):
    """A discone designed from 100 MHz was swept 0.7-1.3 GHz, the 1 GHz default."""
    model = build(registry["discone"].synthesize(f_low=100e6))
    assert model.band_hz == pytest.approx((80e6, 400e6))
    assert 'Solver.FrequencyRange "0.08", "0.4"' in cst.render(model)


def test_probes_stand_clear_of_the_ground_with_the_port_in_the_gap(registry):
    """The DRA probes stood on the unbroken ground (a short) with their ports hanging
    below it in empty space."""
    for key in ("rectangular_dra", "cylindrical_dra", "hemispherical_dra"):
        a = registry[key]
        model = build(a.synthesize(**a.spec.known_cases[0].given))
        probe = next(s for s in model.solids if s.name == "probe")
        port = model.ports[0]
        assert probe.span[0] > 0, key
        assert port.start[2] == 0.0 and port.end[2] == pytest.approx(probe.span[0]), key


def test_waveguides_have_walls(registry):
    """The guides were a vacuum brick (and one slotted wall) inside open boundaries:
    nothing guided anything."""
    for key in ("open_ended_waveguide", "waveguide_slot_array_resonant", "leaky_wave_line_source"):
        a = registry[key]
        names = {s.name for s in build(a.synthesize(**a.spec.known_cases[0].given)).solids}
        assert {"broad_wall", "broad_wall_lower", "narrow_wall_plus", "narrow_wall_minus"} <= names, key
    a = registry["waveguide_slot_array_resonant"]
    assert "end_short" in {s.name for s in build(a.synthesize(**a.spec.known_cases[0].given)).solids}


def test_the_leaky_wave_slit_is_off_the_centreline(registry):
    """TE10's broad-wall current on the centreline is purely longitudinal: a slit
    there cuts none of it and leaks nothing."""
    a = registry["leaky_wave_line_source"]
    cut = next(s for s in build(a.synthesize(**a.spec.known_cases[0].given)).solids if s.name == "slit_cut")
    assert abs(0.5 * (cut.x[0] + cut.x[1])) > (cut.x[1] - cut.x[0])


def test_designed_features_are_built(registry):
    """The inset notch was drawn but never cut; the circular patch's probe ignored
    rho_frac; a three-conductor folded dipole had two; a dielectric spacer in a
    stacked patch was air."""
    from otahub.export.base import Subtract
    inset = build(registry["rectangular_patch_inset"].synthesize(f0=2.4e9, eps_r=4.4, h=0.0016, Z_target=50))
    assert Subtract("patch", ("inset_notch",)) in inset.operations
    for rho in (0.3, 0.7):
        d = registry["circular_patch"].synthesize(f0=10e9, eps_r=2.2, h=0.001588, rho_frac=rho)
        assert build(d).ports[0].start[0] == pytest.approx(d.get("probe_radius_m"))
    three = build(registry["folded_dipole"].synthesize(f0=300e6, N=3, aw=0.001))
    assert sum(1 for s in three.solids if s.name.startswith(("fed_upper", "parasitic"))) == 3
    stack = build(registry["stacked_patch"].synthesize(f0=2.4e9, eps_r=2.2, h=0.0016, eps_r2=2.2))
    upper = next(s for s in stack.solids if s.name == "upper_substrate")
    assert upper.material.startswith("eps_r=2.2") and upper.z[1] - upper.z[0] == pytest.approx(stack.parameters["h2"])


def test_hfss_draws_a_vertical_sheet_upright(registry):
    """The shorting walls are normal to Y; drawn as XY rectangles they had no height."""
    a = registry["quarter_wave_shorted_patch"]
    text = hfss.render(build(a.synthesize(**a.spec.known_cases[0].given)))
    block = text[text.index('"Name:=", "shorting_wall"') - 400:text.index('"Name:=", "shorting_wall"')]
    assert '"WhichAxis:=", "Y"' in block

"""3-D solid models and STL: every archetype gets a closed, consistent solid, its
key dimensions are the design's, and the construction options do what they say.
"""
from __future__ import annotations

import math
import struct
from collections import Counter

import numpy as np
import pytest

mesh = pytest.importorskip("otahub.export.mesh")
if not mesh.available():                                  # pragma: no cover
    pytest.skip("manifold3d is not installed", allow_module_level=True)

from otahub.cli.main import main as cli  # noqa: E402

C0 = 2.99792458e8


def _built(registry, key, options=None):
    """The first known case whose design gives a complete solid."""
    a = registry[key]
    for case in a.spec.known_cases:
        d = a.synthesize(**case.given)
        s = mesh.solid(d, options)
        if s.built:
            return d, s
    raise AssertionError(f"no known case of {key} builds a solid")


def _read_stl(path):
    data = path.read_bytes()
    n = struct.unpack("<I", data[80:84])[0]
    rec = np.frombuffer(data[84:], dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")], count=n)
    return rec["v"].astype(float)


def _closed(tris) -> bool:
    """Every edge shared by exactly two triangles, once in each direction: closed and
    consistently oriented (a duplicated shell fails, which counting alone did not catch)."""
    key = lambda p: tuple(np.round(p, 5))
    edges = Counter()
    for t in tris:
        for i in range(3):
            edges[(key(t[i]), key(t[(i + 1) % 3]))] += 1
    return all(c == 1 and edges.get((b, a), 0) == 1 for (a, b), c in edges.items())


def test_every_archetype_has_a_closed_solid(registry):
    for a in registry:
        d, s = _built(registry, a.key)
        mats = s.materials()
        assert mats, a.key
        for m, solid in mats.items():
            assert solid.status() == mesh.m3.Error.NoError and solid.volume() > 0, (a.key, m)
            assert np.all(np.isfinite(mesh._triangles(solid))), (a.key, m)


@pytest.mark.parametrize("key", ["pyramidal_horn", "rectangular_patch_inset", "axial_mode_helix", "vivaldi_tsa",
                                 "cassegrain", "archimedean_spiral", "hemispherical_dra"])
def test_the_stl_files_are_closed_meshes_in_millimetres(registry, key, tmp_path):
    d, s = _built(registry, key)
    files = mesh.write_stl(s, tmp_path / f"{key}.stl")
    mats = s.materials()
    assert len(files) == (1 if len(mats) == 1 else len(mats) + 1)
    lo, hi = s.bounds()
    for f in files[1:] or files:
        tris = _read_stl(f)
        assert len(tris) > 0 and _closed(tris), f.name
    allt = _read_stl(files[0])
    assert allt.reshape(-1, 3).min(0) == pytest.approx(lo * 1000, abs=1e-3 * max(hi - lo) * 1000)
    assert allt.reshape(-1, 3).max(0) == pytest.approx(hi * 1000, abs=1e-3 * max(hi - lo) * 1000)


def test_the_horn_is_the_design_and_its_wall_option_works(registry):
    d, s = _built(registry, "pyramidal_horn")
    wall = mesh.option_values(d)["wall"]
    lo, hi = s.bounds()
    assert hi[0] - lo[0] == pytest.approx(d.get("a1") + 2 * wall, rel=1e-6)
    assert hi[1] - lo[1] == pytest.approx(d.get("b1") + 2 * wall, rel=1e-6)
    thick = mesh.solid(d, mesh.Options({"wall": 3e-3}))
    lo2, hi2 = thick.bounds()
    assert hi2[0] - lo2[0] == pytest.approx(d.get("a1") + 6e-3, rel=1e-6)
    longer = mesh.solid(d, mesh.Options({"guide_length": 0.1}))
    assert longer.bounds()[0][2] == pytest.approx(-0.1, rel=1e-6)


def test_copper_thickness_and_board_margin_are_the_users(registry):
    a = registry["rectangular_patch_inset"]
    d = a.synthesize(f0=2.4e9, eps_r=3.66, h=0.508e-3)
    s = mesh.solid(d, mesh.Options({"copper": 18e-6, "margin": 0.01}))
    lo, hi = s.bounds()
    assert hi[0] - lo[0] == pytest.approx(d.get("W") + 0.02, rel=1e-6)
    copper = s.materials()["PEC"]
    plain = mesh.solid(d)
    assert copper.volume() < plain.materials()["PEC"].volume() * 0.6     # 18 um of copper against 35
    assert not any("Substrate extended" in n for n in s.notes)


@pytest.mark.parametrize("key,check", [
    ("prime_focus_parabolic", lambda d, lo, hi: hi[0] - lo[0] == pytest.approx(d.get("D"), rel=1e-3)),
    ("yagi_uda", lambda d, lo, hi: hi[0] - lo[0] == pytest.approx(d.get("boom_length"), rel=0.02)),
    ("lpda", lambda d, lo, hi: hi[1] - lo[1] == pytest.approx(d.get("L_max"), rel=0.02)),
    ("axial_mode_helix", lambda d, lo, hi: hi[0] - lo[0] == pytest.approx(d.get("D_gnd"), rel=1e-3)),
    ("half_wave_dipole", lambda d, lo, hi: hi[2] - lo[2] == pytest.approx(d.get("L"), rel=0.01)),
    ("hyperbolic_dielectric_lens", lambda d, lo, hi: hi[0] - lo[0] == pytest.approx(d.get("D_ap"), rel=0.01)),
])
def test_key_dimensions_are_the_designs(registry, key, check):
    d, s = _built(registry, key)
    lo, hi = s.bounds()
    assert check(d, lo, hi)


def test_a_design_without_a_dimension_makes_no_solid(registry):
    d = registry["pifa"].synthesize(f0=2.4e9, h=0.006, W=0.02, Ws=0.001)      # L is NaN here
    s = mesh.solid(d)
    assert not s.built and any("NaN" in n or "unavailable" in n for n in s.notes)


def test_every_option_says_what_it_does_to_the_predictions(registry):
    for a in registry:
        for o in mesh.construction_options(a.key):
            assert o.effect and o.label and o.unit in ("m", "-", "deg"), (a.key, o.name)


def test_the_preview_is_coarser_but_the_same_size(registry):
    d, full = _built(registry, "offset_parabolic")
    with mesh.preview():
        quick = mesh.solid(d)
    assert sum(s.num_tri() for s in quick.materials().values()) < sum(s.num_tri() for s in full.materials().values())
    assert np.allclose(quick.bounds()[1] - quick.bounds()[0], full.bounds()[1] - full.bounds()[0], rtol=0.02)


def test_the_cli_writes_stl_and_lists_options(tmp_path, capsys):
    out = tmp_path / "patch.stl"
    assert cli(["export", "rectangular_patch_inset", "--f0", "2.4GHz", "--set", "eps_r=4.4", "--set", "h=1.6mm",
                "--format", "stl", "--opt", "copper=35um", "-o", str(out)]) == 0
    assert out.exists() and (tmp_path / "patch_conductor.stl").exists()
    assert cli(["export", "pyramidal_horn", "--f0", "10GHz", "--set", "G_target=20", "--options"]) == 0
    text = capsys.readouterr().out
    assert "wall" in text and "guide_length" in text and "geometry only" in text
    assert cli(["export", "pyramidal_horn", "--f0", "10GHz", "--set", "G_target=20", "--format", "stl",
                "--opt", "nonsense=1", "-o", str(tmp_path / "h.stl")]) == 2


@pytest.fixture(scope="module")
def qapp():
    widgets = pytest.importorskip("PySide6.QtWidgets")
    return widgets.QApplication.instance() or widgets.QApplication([])


def test_the_gui_shows_the_model_and_saves_it(qapp, registry, tmp_path):
    from otahub.gui.model3d import SolidPanel
    panel = SolidPanel()
    d, _ = _built(registry, "conical_horn")
    panel.set_design(d)
    panel.refresh()
    assert "Conical horn" in panel.summary.text()
    assert set(panel._fields) == {o.name for o in mesh.construction_options("conical_horn")}
    panel._fields["wall"].setText("2 mm")
    files = panel.save_stl(str(tmp_path / "horn.stl"))
    assert files and files[0].exists()
    lo, hi = mesh.solid(d, mesh.Options({"wall": 2e-3})).bounds()
    tris = _read_stl(files[0]).reshape(-1, 3)
    assert tris.max(0) - tris.min(0) == pytest.approx((hi - lo) * 1000, rel=1e-4)


def test_a_duplicated_shell_is_not_closed():
    tet = np.array([[[0, 0, 0], [0, 1, 0], [1, 0, 0]], [[0, 0, 0], [1, 0, 0], [0, 0, 1]],
                    [[0, 0, 0], [0, 0, 1], [0, 1, 0]], [[1, 0, 0], [0, 1, 0], [0, 0, 1]]], float)
    assert _closed(tet) and not _closed(np.concatenate([tet, tet]))


# ---------------------------------------------------------------- the review's findings

@pytest.mark.parametrize("key", ["triangular_patch", "truncated_corner_cp_patch"])
def test_a_patch_probe_does_not_short_the_patch_to_ground(registry, key):
    d = registry[key].synthesize(f0=2.4e9, eps_r=2.2, h=0.0016)
    pieces = mesh.solid(d).materials()["PEC"].decompose()
    assert len(pieces) == 2                     # ground; patch + probe through its clearance hole


def test_the_horn_hollow_is_the_designed_aperture_and_throat(registry):
    d = registry["pyramidal_horn"].synthesize(f0=10e9, G_target=20)
    horn = mesh.solid(d).bodies[0].solid
    for z, a, b in ((d.get("p_len") - 1e-7, d.get("a1"), d.get("b1")), (1e-7, d.get("a_wg"), d.get("b_wg"))):
        inner = min((np.asarray(q) for q in horn.slice(z).to_polygons()), key=lambda q: np.ptp(q[:, 0]))
        # the mesh keeps coordinates in single precision, and the slice sits 0.1 um into the flare
        assert np.ptp(inner[:, 0]) == pytest.approx(a, rel=1e-5) and np.ptp(inner[:, 1]) == pytest.approx(b, rel=1e-5)


def test_a_finite_thin_slot_keeps_its_width(registry):
    d = registry["half_wave_slot"].synthesize(f0=100e9, w_over_L=0.05)
    plate_ = mesh.solid(d).materials()["PEC"]
    polys = plate_.slice(-1e-7).to_polygons()
    hole = min((np.asarray(q) for q in polys), key=lambda q: np.ptp(q[:, 0]) * np.ptp(q[:, 1]))
    assert min(np.ptp(hole[:, 0]), np.ptp(hole[:, 1])) == pytest.approx(d.get("w"), rel=1e-4)    # single precision


def test_the_margin_covers_a_top_hat(registry):
    a = registry["top_loaded_monopole"]
    d = a.synthesize(**a.spec.known_cases[0].given)
    s = mesh.solid(d, mesh.Options({"margin": 0.01}))
    ground = next(b for b in s.bodies if b.name in ("ground", "ground_plane", "gnd")).solid.bounding_box()
    others = [b.solid.bounding_box() for b in s.bodies if b.name not in ("ground", "ground_plane", "gnd")]
    assert ground[3] >= max(o[3] for o in others) + 0.01 - 1e-9


def test_the_vivaldi_opens_at_the_designs_rate(registry):
    v = [mesh.solid(registry["vivaldi_tsa"].synthesize(f_low=1e9, Lax=0.3, R_open=R)).materials()["PEC"].volume()
         for R in (10, 20)]
    assert v[1] > v[0] * 1.05                   # a faster-opening slot leaves more copper near the throat


def test_an_explicit_feed_size_wins_over_the_blockage_default(registry):
    d = registry["prime_focus_parabolic"].synthesize(f0=10e9, D=1, d_blockage=0.1)
    widths = []
    for fd in (0.03, 0.2):
        feed = next(b for b in mesh.solid(d, mesh.Options({"feed_diameter": fd})).bodies if b.name == "feed")
        bb = feed.solid.bounding_box()
        widths.append(bb[3] - bb[0])
    assert widths[1] - widths[0] == pytest.approx(0.17, rel=1e-3)


def test_the_ferrite_coil_keeps_its_length_and_centre(registry):
    d = registry["ferrite_rod_loop"].synthesize(f0=1e6, N=60, l_coil=0.02, b=0.0002, l_rod=0.1, d_rod=0.01, mu_i=125)
    coils = [b.solid.bounding_box() for b in mesh.solid(d).bodies if b.name.startswith("coil")]
    assert len(coils) == 2                      # 60 turns of 0.4 mm wire need two layers on 20 mm
    for bb in coils:
        assert bb[0] == pytest.approx(-bb[3], abs=1e-6) and bb[3] - bb[0] <= 0.02 + 2 * 0.0002 + 1e-6


def test_the_long_wire_sits_at_the_designs_height(registry):
    d = registry["long_wire_travelling"].synthesize(f0=100e6, L_over_lambda=4, height=1.0, aw=0.001)
    assert mesh.option_values(d)["height"] == pytest.approx(1.0)
    assert mesh.solid(d).bounds()[1][2] == pytest.approx(1.0 + 0.001, rel=1e-9)


def test_the_pifa_strip_is_at_a_corner(registry):
    from otahub.export import build
    d = registry["pifa"].synthesize(f0=1e9, h=0.012, W=0.036, Ws=0.006)
    m = build(d)
    wall = next(s for s in m.solids if s.name == "shorting_wall")
    assert wall.x == pytest.approx((-0.018, -0.012))
    assert m.ports[0].start[0] == pytest.approx(-0.015)


@pytest.mark.parametrize("key, opt, value, says", [
    ("corrugated_conical_horn", "pitch", 0.0, "positive"),
    ("rectangular_patch_inset", "copper", 0.0, "positive"),
    ("rectangular_patch_inset", "copper", -1.0, "positive"),
    ("rectangular_patch_inset", "margin", -0.03, "zero or more"),
    ("corrugated_conical_horn", "tooth", 1.5, "between 0 and 1"),
    ("prime_focus_parabolic", "struts", 2.5, "count"),
    ("prime_focus_parabolic", "struts", float("inf"), "finite"),
])
def test_an_unusable_option_is_refused_with_a_reason(registry, key, opt, value, says):
    """A zero pitch looped forever (the GUI froze); copper of zero or below dropped every
    conductor from the STL while the CLI reported success; a negative margin cut the board
    under the patch. Each is now a note, and no model."""
    assert any(says in p for p in mesh.option_problems(key, {opt: value}))
    case = registry[key].spec.known_cases[0]
    s = mesh.solid(registry[key].synthesize(**case.given), mesh.Options({opt: value}))
    assert not s.built and says in s.notes[0]


def test_a_pitch_too_fine_to_build_stops_at_once(registry):
    d = registry["corrugated_conical_horn"].synthesize(f0=10e9, L=0.3, flare=0.1)
    s = mesh.solid(d, mesh.Options({"pitch": 3.7e-9}))
    assert not s.built and "at most 2000" in s.notes[0]


def test_option_defaults_survive_an_incomplete_design(registry):
    """`--options` on an E-plane horn without its guide width raised KeyError: 'a_wg'
    (its b_wg option's default reads it), and in the GUI left the previous horn on screen."""
    d = registry["e_plane_sectoral_horn"].synthesize(f0=10e9)
    values = mesh.option_values(d)
    assert math.isnan(values["b_wg"]) and values["wall"] > 0


def test_a_loop_without_its_wire_radius_names_it(registry):
    for key, f0 in (("halo_loop", 50.1e6), ("one_wavelength_circular_loop", 3e8), ("quad_loop_square", 14.2e6)):
        s = mesh.solid(registry[key].synthesize(f0=f0))
        assert not s.built and "waiting for b" in s.notes[0] and "('" not in s.notes[0]

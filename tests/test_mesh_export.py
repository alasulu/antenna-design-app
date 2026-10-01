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
    """Every edge shared by exactly two triangles, in opposite directions."""
    key = lambda p: tuple(np.round(p, 5))
    edges = Counter()
    for t in tris:
        for i in range(3):
            edges[(key(t[i]), key(t[(i + 1) % 3]))] += 1
    return all(edges.get((b, a), 0) == c for (a, b), c in edges.items())


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

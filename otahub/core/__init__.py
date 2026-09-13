"""Core engine: GUI-free, depends only on numpy/scipy.

Import nothing from a GUI toolkit here — the CLI, the test suite and any
future headless batch runs all rely on this staying importable on a bare
scientific Python install.
"""
from .archetype import Archetype, DesignResult, SynthesisError
from .constants import C0, EPS0, ETA0, MU0, skin_depth, surface_resistance, wavelength, wavenumber
from .expr import ExprError, evaluate
from .registry import Registry, default_registry
from .spec import ArchetypeSpec, FamilySpec, KnownCase

__all__ = [
    "Archetype", "DesignResult", "SynthesisError",
    "ArchetypeSpec", "FamilySpec", "KnownCase",
    "Registry", "default_registry",
    "evaluate", "ExprError",
    "C0", "EPS0", "MU0", "ETA0",
    "wavelength", "wavenumber", "skin_depth", "surface_resistance",
]

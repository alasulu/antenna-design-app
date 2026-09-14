"""Array layout, excitation and pattern synthesis."""
from .factor import (array_factor, array_pattern, beam_broadening_factor,
                     broadside_directivity, grating_lobe_free_spacing,
                     has_grating_lobe, summarise)
from .planar import (grating_lobe_free_spacing_planar, lattice_element_saving,
                     planar_beam_cut,
                     planar_array_factor, planar_directivity, planar_pattern_cut,
                     planar_summarise, rectangular_lattice, separable_weights,
                     steering_phase, triangular_lattice)
from .tapers import (TAPERS, binomial, dolph_chebyshev, raised_cosine,
                     taper_efficiency, taylor_nbar, uniform)

__all__ = ["uniform", "binomial", "dolph_chebyshev", "taylor_nbar",
           "raised_cosine", "taper_efficiency", "TAPERS",
           "array_factor", "array_pattern", "broadside_directivity",
           "grating_lobe_free_spacing", "has_grating_lobe",
           "beam_broadening_factor", "summarise",
           "rectangular_lattice", "triangular_lattice", "separable_weights",
           "steering_phase", "planar_array_factor", "planar_directivity",
           "planar_pattern_cut", "planar_beam_cut", "grating_lobe_free_spacing_planar",
           "lattice_element_saving", "planar_summarise"]

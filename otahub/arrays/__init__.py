"""Array layout, excitation and pattern synthesis."""
from .elements import (ElementPattern, cosine, custom, element_directivity, isotropic,
                       power_kernel, short_dipole)
from .layouts import (circular_taylor, clip_to_circle, concentric_rings, equal_area_radius, ring,
                      taylor_circular_distribution, taylor_circular_pattern, thin,
                      thinned_expected_directivity, thinned_expected_power)
from .subarrays import directivity_toward, rect_subarray_groups, subarray_centres, subarray_steering
from .factor import (array_factor, array_pattern, beam_broadening_factor,
                     broadside_directivity, directivity, grating_lobe_free_spacing,
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
           "array_factor", "array_pattern", "broadside_directivity", "directivity",
           "grating_lobe_free_spacing", "has_grating_lobe",
           "beam_broadening_factor", "summarise",
           "rectangular_lattice", "triangular_lattice", "separable_weights",
           "steering_phase", "planar_array_factor", "planar_directivity",
           "planar_pattern_cut", "planar_beam_cut", "grating_lobe_free_spacing_planar",
           "lattice_element_saving", "planar_summarise",
           "ElementPattern", "isotropic", "cosine", "short_dipole", "custom",
           "power_kernel", "element_directivity",
           "ring", "concentric_rings", "clip_to_circle", "taylor_circular_pattern",
           "taylor_circular_distribution", "equal_area_radius", "circular_taylor", "thin",
           "thinned_expected_power", "thinned_expected_directivity",
           "rect_subarray_groups", "subarray_centres", "subarray_steering", "directivity_toward"]

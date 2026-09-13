"""Array layout, excitation and pattern synthesis."""
from .factor import (array_factor, array_pattern, beam_broadening_factor,
                     broadside_directivity, grating_lobe_free_spacing,
                     has_grating_lobe, summarise)
from .tapers import (TAPERS, binomial, dolph_chebyshev, raised_cosine,
                     taper_efficiency, taylor_nbar, uniform)

__all__ = ["uniform", "binomial", "dolph_chebyshev", "taylor_nbar",
           "raised_cosine", "taper_efficiency", "TAPERS",
           "array_factor", "array_pattern", "broadside_directivity",
           "grating_lobe_free_spacing", "has_grating_lobe",
           "beam_broadening_factor", "summarise"]

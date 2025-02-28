module Erbium


# local functions
export normalize_spectrum
export power_spectrum
export mean_wavelength
export bandwidth
export dbm2mw

export optical_fiber
export generate_gaussian_spectrum
export dIdz
export dIpdz
export dIdz_opt
export dIpdz_opt


# bandiwdth() calculation types to allow multiple dispatch
abstract type bwMethod end
struct FWHM <: bwMethod end
struct Weighted <: bwMethod end

export bwMethod, FWHM, Weighted


# 3rd party package functions
export linear_interpolation
export trapz

using DelimitedFiles
using Interpolations
using Trapz


include("rate_equations.jl")
include("optical_fiber.jl")
include("spectral_analysis.jl")
include("physics_constants.jl")




end # module Erbium

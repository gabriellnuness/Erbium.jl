module Erbium


# local functions
export normalize_spectrum
export power_spectrum
export mean_wavelength
export bandwidth
export dbm2mw
export test_calc_spec_parameters
export test_integrator_single_pass

export optical_fiber
export generate_gaussian_spectrum
export dPdz
export dIpdz
export dIdz_opt
export dIpdz_opt
export trapz
export pump_power_980
export pump_power_1480

# bandiwdth() calculation types to allow multiple dispatch
abstract type bwMethod end
struct FWHM <: bwMethod end
struct Weighted <: bwMethod end

export bwMethod, FWHM, Weighted


# 3rd party package functions
export linear_interpolation

using DelimitedFiles
using Interpolations

include("rate_equations.jl")
include("optical_fiber.jl")
include("spectral_analysis.jl")




end # module Erbium

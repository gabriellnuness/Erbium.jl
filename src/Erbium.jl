module Erbium

# local functions
export optical_fiber
export generate_gaussian_spectrum
export dIdz
export dIpdz

# 3rd party package functions
export linear_interpolation
export trapz

using DelimitedFiles
using Interpolations
using Trapz


include("rate_equations.jl")
include("optical_fiber.jl")




end # module Erbium

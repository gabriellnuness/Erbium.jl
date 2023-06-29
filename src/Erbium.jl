module Erbium

export dJpump_dz
export dJ_dz
export optical_fiber


using DelimitedFiles
using Interpolations
using Trapz


include("rate_equations.jl")
include("optical_fiber.jl")




end # module Erbium

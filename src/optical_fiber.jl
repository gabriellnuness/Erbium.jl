"""
Set constants according to chosen optical fiber.
(β_abs, β_emis, τ21, τ3, diameter, NA, total_population, η) = optical_fiber(fiber)

#parameters:
* `fiber`: optical fiber model `String` can be `m5` or `m12`.

#returns:
* `β_abs`, `β_emis`:        Absorption and emission spectrum provided by manufacturer or from experiment.
* `τ21`, `τ3`:              Life time of erbium levels.
* `diameter`:               Optical mode diameter. 
* `NA`:                     Numerical aperture.
* `total_population`:       Ions population in 1/m³  (Guess!)
* `η`:                      Fluorescence efficience
"""
function optical_fiber(fiber::String)

    # Optical fiber diameter
    if fiber == "m5"
        
        # β = σ⋅NT
        β_abs = readdlm("data\\M5_abs.txt",',')
        β_emis = readdlm("data\\M5_emis.txt",',')
        diameter = ((5.7+6.6) / 2 )*1e-6
        total_population = 2.84e24
        NA = 0.24
        τ21 = 10e-3
        τ3 = 1e-6
        η = 0.1

        return (β_abs, β_emis, τ21, τ3, diameter, NA, total_population, η)

    elseif fiber == "m12"
    
        β_abs = readdlm("data\\M12_abs.txt",',')
        β_emis = readdlm("data\\M12_emis.txt",',')
        diameter = ((5.5+6.3) / 2)*1e-6
        total_population = 2.84e24
        τ21 = 1e-6
        τ3 = 10e-3
        NA = 0.24
        η = 0.1
    
        return (β_abs, β_emis, τ21, τ3, diameter, NA, total_population, η)
        
    
    else
        throw(error("Fiber model not found."))
    end


end


"""
Generate gaussian spectrum in order to represent the pump
"""
function generate_gaussian_spectrum(λ, λc, Δλ, I_peak)
    dλ = λ[2]-λ[1]
    I_pump = @.  √(2/π)*I_peak/dλ*exp(-(((λ-λc)/Δλ)^2))

    return I_pump
end

"""
Distribution of intensity inside optical fiber

Calculated so that the optical fiber diameter is to 3ω₀, 
where ω₀ is the Gaussian linewidth of the distribution.

### Arguments
- `λ`: Wavelength
"""
# function Γ(λ::Float32)
#     w = 0.66*λ # ?
#     a = 1 # ?
#     1-exp(-2*(a/w)^2) # should it be 0.66?
# end


"""
Photon intensity

g(λ) represents the σₑ, spectrum of spontaneous emission.
"""
# function p(λ::Float32)
#     g = λ*1 # import emission spectrum here
#     NA = 0.225  # numerical apperture optical fiber
#     h*v*g*η/2*(1-(1-NA^2)^(1/2))
# end

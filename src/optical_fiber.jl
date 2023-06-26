using DelimitedFiles

"""
Set constants according to chosen optical fiber.

### Arguments
`fiber::String`: optical fiber model can be `m5` or `m12`.
"""
function optical_fiber(fiber::String)

    # Optical fiber diameter
    if fiber == "m5"
        
        # β = σ⋅NT
        β_abs = readdlm("data\\M5_abs.txt",',')
        β_emis = readdlm("data\\M5_emis.txt",',')
        
        diameter = (5.7+6.6) / 2
        # sigma_p_a = 
        # sigma_p_e = 
        # sigma_a = 
        # sigma_e = 
        # N2 = 
        # NT = 
        τ21 = 1e-6
        τ3 = 10e-3

        (β_abs, β_emis, τ21, τ3, diameter)

    elseif fiber == "m12"
    
        diameter = (5.5+6.3) / 2
        # sigma_p_a = 
        # sigma_p_e = 
        # N2 = 
        # NT = 
        # sigma_a = 
        # sigma_e = 
        τ21 = 1e-6
        τ3 = 10e-3
        β_abs = readdlm("data\\M12_abs.txt",',')
        β_emis = readdlm("data\\M12_emis.txt",',')
    
        (β_abs, β_emis, τ21, τ3, diameter)
    end

end


"""
Distribution of intensity inside optical fiber

Calculated so that the optical fiber diameter is to 3ω₀, 
where ω₀ is the Gaussian linewidth of the distribution.

### Arguments
- `λ`: Wavelength
"""
function Γ(λ)
    w = 0.66*λ # ?
    a = 1 # ?
    1-exp(-2*(a/w)^2) # should it be 0.66?
end


"""
Photon intensity

g(λ) represents the σₑ, spectrum of spontaneous emission.
"""
function p(λ)
    g = λ*1 # import emission spectrum here
    NA = 0.225  # numerical apperture optical fiber
    h*v*g*η/2*(1-(1-NA^2)^(1/2))
end

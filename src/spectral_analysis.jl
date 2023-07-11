"""
Δλ = bandwidth(λ,y)

Function to calculate the spectral bandwidth of a given spectrum given in dBm.
This function takes into consideration the double-peaked shape of the
erbium broadband emission. The 3-dB FWHM won't work in this case, so the
appropriate definition is weighted by the square of the power.

# arguments:
- `[λ]` and `[y]`: vectors of λ and y.

# returns:
- `Δλ`: bandwidth value.

# ref:
Falquier, D. G.; 2000. Stanford.
"""
function bandwidth(λ,y)
    
    dλ = diff(λ)
    Δλ = @. sum(y[2:end]*dλ)^2 / sum((y[2:end]^2)*dλ)

    return Δλ
    
end








"""
λ_mean = mean_wavelength(λ,y)

Calculates the mean wavelength of a spectrum.

# arguments:
- `[λ]` and `[y]`: vectors of λ and y.

# returns:
- `λ_mean`: mean wavelength value.

# ref:
Falquier, D. G.; 2000. Stanford.
"""
function mean_wavelength(λ,y)
    
    dλ = diff(λ)
    λ_mean = @. sum(y[2:end]*λ[2:end]*dλ) / sum(y[2:end]*dλ)

    return λ_mean
end








"""
power = spectrum_power(λ,y)

Calculates the integrated power of an spectrum.

# arguments:
- `[λ]` and `[y]`: vectors of λ and y.

# returns:
- `λ_mean`: mean wavelength value.

# warnings:
- Normalize the optical spectrum power to 1 nm resolution 
before using as the function arguments 

"""
function spectrum_power(λ,y)

    dλ = diff(λ); 
    power = @. sum(y[2:end]*dλ);

    return power
end






"""
normalized_power = normalize_spectrum(λ, y)

Normalization of spectrum from Advantest Optical Spectrum Analizer.
The OSA provides the ∫P⋅dλ, where dλ is the OSA resolution for that measurement.

# arguments:
- power [dBm]
- λ     [m]

# returns:
- normalized spectrum [dBm]

"""           
function normalize_spectrum(λ, y)

    dλ = diff(λ)
    normalized_power = @. y[2:end] - 10*log10(dλ*1e9);

    return normalized_power
end
#TODO: Implement tests for these functions
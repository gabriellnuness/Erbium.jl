"""
Δλ = bandwidth(λ, y, ::bwMethod, db_check)

Function to calculate the spectral bandwidth of a given spectrum given in linear units.
This function takes into consideration the double-peaked shape of the
erbium broadband emission. The FWHM won't work in this case, so the
appropriate definition is weighted by the square of the power.

# arguments:
- `[λ]` and `[y]`:  vectors of λ and y [linear].
- `::bwMethod`:     can be `FWHM` or `Weighted`.
- `db_check`:       calculates with curve in decibel if equals "dB".

# returns:
- `Δλ`: bandwidth value.

# ref:
Falquier, D. G.; 2000. Stanford.
"""
function bandwidth(λ, y, ::Type{Weighted})
    
    dλ = diff(λ)
    Δλ = sum(y[2:end].*dλ)^2 / sum((y[2:end].^2).*dλ)

    return Δλ
end
"""
Δλ = bandwidth(λ, y, method)

Function to calculate the spectral bandwidth of a given spectrum using 
the full width at half maximum method.

# arguments:
- `[λ]` and `[y]`: vectors of λ and y [linear].

# returns:
- `Δλ`: bandwidth value.
"""
function bandwidth(λ, y, ::Type{FWHM}, db_check=nothing)
   
    ymax = maximum(y)
    if db_check == "dB"
        y_half = ymax-3 
    elseif isnothing(db_check)
        y_half = ymax/2
    end

    ind1 = findfirst(>=(y_half), y)
    ind2 = findlast(>=(y_half), y)
    λ₁ = λ[ind1]
    λ₂ = λ[ind2]

    # check if the points are an exact match
    # it depends on the resolution of the input data
    # then, make a linear fit to find exact λ representing middle point
    if λ₁ ≉  y_half
        ind12 = ind1 - 1
        X = [λ[ind1], λ[ind12]] # 1x2
        Y = [y[ind1], y[ind12]] # 1x2
        A = [[1,1] X]
        sol = A \ Y    # 2x1 
        m = sol[2]
        b = sol[1]
        
        λ₁ = (y_half-b)/m
    end
    if λ₂ ≉  y_half
        ind22 = ind2 + 1
        X = [λ[ind2], λ[ind22]]
        Y = [y[ind2], y[ind22]]
        A = [[1,1] X]
        sol = A \ Y 
        m = sol[2]
        b = sol[1]
        
        λ₂ = (y_half-b)/m
    end


    Δλ = λ₂ - λ₁
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
    λ_mean = sum(y[2:end].*λ[2:end].*dλ) / sum(y[2:end].*dλ)

    return λ_mean
end




"""
Converts dBm to mW
"""
dbm2mw(dbm) = 10 .^(dbm/10)




"""
power = power_spectrum(λ,y)

Calculates the integrated power of an spectrum.

# arguments:
- `[λ]`: vector of wavelentgh.
- `[y]`: vectors of optical power in Watts.

# returns:
- `λ_mean`: mean wavelength value.

# warnings:
- Normalize the optical spectrum power 
before using as the function argument.

"""
function power_spectrum(λ,y)

    dλ = diff(λ)
    power = sum(y[2:end].*dλ)

    return power
end







"""
normalized_power = normalize_spectrum(λ, y)

Normalization of spectrum from Advantest Optical Spectrum Analizer.
The OSA provides the ∫P⋅dλ, where dλ is the OSA resolution for that measurement.

# arguments:
- `power`:      Vector of power in Watts.
- `λ`:          Vector of wavelength in meters.
- `dλ`:         Float value for the resolution. For the cases in which the OSA provides a fixed resolution that can't be obtained from derivating `λ`.

# returns:
- normalized spectrum [W]

"""           
function normalize_spectrum(λ, y)

    dλ = diff(λ)
    normalized_power = y[2:end] ./ dλ

    return (λ_norm=λ[2:end], power_norm=normalized_power)
end
function normalize_spectrum(λ, y, dλ)

    normalized_power = y/dλ

    return (λ_norm=λ, power_norm=normalized_power)
end

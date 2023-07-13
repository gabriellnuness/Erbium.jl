using Erbium
using PyPlot
using DelimitedFiles


power_ref = 4.74e-3  # reference power measured with power meter
input = readdlm("data\\example_sld_4_74mw.lvm", ',')
power_dbm = input[:,1]
λ = input[:,2]


# Convert dBm to Watts
power = dbm2mw(power_dbm)*1e-3
# Normalize the spectrum remove dλ dependence
λ_norm, power_norm = normalize_spectrum(λ, power)

# Perform all calculations with normalized optical power in SI
total_power = power_spectrum(λ_norm, power_norm)
λ_mean = mean_wavelength(λ_norm, power_norm)

Δλ = bandwidth(λ_norm, power_norm, Weighted)
Δλ = bandwidth(λ_norm, power_norm, FWHM)


# normalized plot compared to the original in W
# figure()
# plot(λ_norm*1e9, power_norm)
# ylabel("Power normalized in Watts  = Pₒₛₐ / dλ")
# xlabel("λ [nm]")
# ax = gca().twinx()
# ax.plot(λ*1e9, power,color="tab:orange")
# ylabel("Power in Watts = Pₒₛₐ")


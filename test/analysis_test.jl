using Erbium
using PyPlot
using DelimitedFiles

##
##
## SLD example
##
##
package_path = dirname(dirname(pathof(Erbium)))
input =readdlm("$(package_path)/data/example_sld_4_74mw.lvm", ',')
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
Δλ = bandwidth(λ, power_dbm, FWHM, "dB")


# normalized plot compared to the original in W
figure()
plot(λ_norm*1e9, power_norm)
ylabel("Power normalized in Watts  = Pₒₛₐ / dλ")
xlabel("λ [nm]")
ax = gca().twinx()
ax.plot(λ*1e9, power,color="tab:orange")
ylabel("Power in Watts = Pₒₛₐ")


##
##
## Erbium example
##
##
input = readdlm("$(package_path)/data/example_erbium_4_44mw.lvm", ',')
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
Δλ = bandwidth(λ, power_dbm, FWHM, "dB")


# normalized plot compared to the original in W
figure()
plot(λ_norm*1e9, power_norm)
ylabel("Power normalized in Watts  = Pₒₛₐ / dλ")
xlabel("λ [nm]")
ax = gca().twinx()
ax.plot(λ*1e9, power,color="tab:orange")
ylabel("Power in Watts = Pₒₛₐ")


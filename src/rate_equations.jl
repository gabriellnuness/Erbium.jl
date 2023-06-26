include("optical_fiber.jl")

# differential equations to simulate broadband amplified erbium spectrum
# ----------- 980 nm / 1480 nm -----------
#   ^    |               |
#   |    |              τ3
#   |    |              v
#   Wpa  Wpe    --- 1530 ≈ 1550 nm -------
#   |    |           ^      |     |
#   |    |          Wa     We     τ21
#   |    v          |      v      v
# ----------------------------------------


# TODO: integrate the pump_absorption in relation to λ as well.
# TODO: define the sigmas as the function of the spectra in \\data folder
pump_absorption(λ) = -Γ(λ) * sigma_p_a(λ) * (NT-N2(z))
pump_stimulated_emission(λ) = Γ(λ) * sigma_p_e(λ) * (NT-N2(z))

beam_absorption(λ) = -Γ(λ) * sigma_a(λ) * (NT-N2(z))
beam_stimulated_emission(λ) = Γ(λ) * sigma_e(λ) * N2(z)
beam_natural_emission(λ) = p(λ) * N2(z)/τ21


# Trapz integration rule
# W = ∫(λ^2*h*v/sigma(λ)*c)dλ


N2(z) = (pump_absorption + beam_absorption)*τ21*NT /
        ((pump_absorption + beam_absorption + beam_stimulated_emission)*τ21 + 1)



"""Pump rate equation"""
dJpump_dz(λ,z) = Jpump(λ,z) * (-pump_absorption + pump_stimulated_emission - gamma(λ))


""" Beam rate equation"""
dJ_dz(λ,z) = J(λ,z) * (+beam_absorption + beam_stimulated_emission - gamma(λ)) + beam_natural_emission


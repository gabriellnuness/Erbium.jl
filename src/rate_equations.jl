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



"""Beam rate equation"""
function dIdz(I,Γ,β_abs,n2,β_emis,fiber_loss,Rho)

    absorption = Γ*β_abs*(1-n2)
    emission = Γ*β_emis*n2
    spontaneous = Rho*n2

    return I*(-absorption + emission - fiber_loss) + spontaneous
end



"""Pump rate equation"""
function dIpdz(I,Γ,β_abs_pump,n2,β_emis_pump,fiber_loss)

    pump_absorption = Γ*β_abs_pump*(1-n2)
    pump_emission =  Γ*β_emis_pump*n2

    return I*(-pump_absorption + pump_emission - fiber_loss)
end



"""
4th order Runge-Kutta integration method to solde the beam rate ODE 
"""
function ∫()
      

end
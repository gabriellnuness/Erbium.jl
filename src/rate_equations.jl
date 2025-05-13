# differential equations to simulate broadband amplified erbium spectrum
# -------------- 980 nm  -------------------------
#   ^    |               |
#   |    |              τ3
#   |    |              v
#   Wpa  Wpe  ----- 1480 ≈ 1530 ≈ 1550 nm --------
#   |    |           ^      |     |
#   |    |          Wa     We     τ21
#   |    v          |      v      v
# ------------------------------------------------



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

# arguments:
* `fun`:            function to be integrated.
* `dx`:             integration step.
* `xk` and `yk`:    initial values for `x[k]` and `y[k]`.

# returns:
* `yout`:           function evaluated in x[k+1].
"""
function rk4_fixedstep(fun, dx, xk, yk)
      f1 = fun(xk,        yk            )
      f2 = fun(xk+dx/2,   yk + (dx/2)*f1)
      f3 = fun(xk+dx/2,   yk + (dx/2)*f2)
      f4 = fun(xk+dx,     yk + dx*f3    )       
      yout = yk + (dx/6)*(f1 + 2*f2 + 2*f3 + f4)
      return yout
end
"""

f1 = dIdz(I_forward[λ,z],           Γ,β_abs[λ],n2[z],β_emis[λ],fiber_loss,Rho)*dz
f2 = dIdz(I_forward[λ,z] + f1/2,    Γ,β_abs[λ],n2[z],β_emis[λ],fiber_loss,Rho)*dz
f3 = dIdz(I_forward[λ,z] + f2/2,    Γ,β_abs[λ],n2[z],β_emis[λ],fiber_loss,Rho)*dz
f4 = dIdz(I_forward[λ,z] + f3,      Γ,β_abs[λ],n2[z],β_emis[λ],fiber_loss,Rho)*dz
I_forward[λ,z+1] = I_forward[λ,z] + (f1+ 2*f2 + 2*f3 + f4)/6


I = rk4fixedstep(dIdz(I,Γ,β_abs,n2,β_emis,fiber_loss,Rho), dz, z, I0)
"""

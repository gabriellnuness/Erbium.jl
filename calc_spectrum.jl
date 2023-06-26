using Erbium

optical_fiber("m5")

J = RK4(dJ_dz, initial_conditions, tspan)
Jpump = RK4(dJpump_dz, initial_conditions, tspan)

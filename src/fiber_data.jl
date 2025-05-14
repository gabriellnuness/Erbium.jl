
# Nz = 300 #Número de Intervalos em z
# Nl = 300# Número de Intervalos de λs
# Nk = 500 # Número máximo de iterações
# Nc = Int(Nl/4)  # Corte no espectro para desconsiderar bomba 
# Nexp=2500
# Tol=0.001 # Tolerância para testes de convergência

# λ0 =1550e-9;
# λ1 = 1450e-9 #nanometros - valor inferior λ absorção/emissão
# λ2 = 1600e-9 #nanometros - valor superior λ absorção/emissão


# #### Dados da fibra
# RL_980 = 0.0 #2.3e-5 #0.00001# # Refletividade no fim da fibra
# RL_1480 = 0.3
# R0_980 = 0.0 #0.9#1e-5
# R0_1480 = 0.0
# T = 273+27  #Temperatura em kelvin
# Ksigma = 0.85 # Datasheet FiberCore para M12 0.9
# Z= 1.29e-9 #  Ver estimativa inicial de Z e Rho no Rho.sm SMath
# ϵ1=0.1
# ϵ2=0.0159
# Rho0 = Z*ϵ1*ϵ2
# gama = 0.01 # perda da fibra em m^-1
# G= 1.0  # Fator de sobreposição entre feixe e distribuição de Er na fibra

# fiber_parameters = (RL_980 = 0.0,
# 					RL_1480 = 0.3,
# 					R0_980 = 0.0,
# 					R0_1480 = 0.0,
# 					T = 273+27,
# 					Ksigma = 0.85,
# 					Z= 1.29e-9,
# 					ϵ1=0.1,
# 					ϵ2=0.0159,
# 					gama = 0.01,
# 					G= 1.0
# 					)
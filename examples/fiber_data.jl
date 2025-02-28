

"""
	begin include("nicolau_Dados-M5-980.jl")
"""

""" Constantes """

h=6.6e-34    	# Constante de Planck
c = 3e8      	# Velocidade da luz
kB = 1.38e-23  	# Constante de Boltzmann
Nz = 300 		# Número de Intervalos em z
Nλ = 300		# Número de Intervalos de λs
max_iter = 500 	# Número máximo de iterações
Nc = Int(Nλ/4)  # Corte no espectro para desconsiderar bomba 
Nexp = 2500
Tol = 0.001 	# Tolerância para testes de convergência (0.1%)

# bandwidth of desired spectrum
λ0 = 1550e-9
λ1 = 1450e-9	# nanometros - valor inferior λ absorção/emissão
λ2 = 1600e-9 	# nanometros - valor superior λ absorção/emissão


""" Dados da fibra """

RL = 2.3e-5 	# 0.00001# # Refletividade no fim da fibra
R0 = 0.9		# 1e-5
T = 273+27  	# Temperatura em kelvin
Ksigma = 0.85 	# Datasheet FiberCore para M12 0.9
Z0 = 1.15
Z= Z0*1.58e-9	#  Ver estimativa inicial de Z e Rho no Rho.sm SMath
ϵ1 = 0.1
ϵ2 = 0.0080 	# 0.0159
Rho0 = Z*ϵ1*ϵ2
gama = 0.0014 	# perda da fibra em m^-1
G = 1.0  		# Fator de sobreposição entre feixe e distribuição de Er na fibra
L = 10 			# metros - comprimento da fibra dopada com érbio ###


""" Dados da bomba """

P0_980  = -0.025+0.65*I_980  # 40e-3 #Potência da bomba em watt ( 980 nm)
if P0_980 <0
	P0_980 = 0 
end
I_1480 = 0 				# Corrente do laser 1480 nm
P0_1480 = 0 			# 0.009+0.153*I_1480 #20e-3 #Potência da bomba em watt (1480 nm)
if P0_1480 < 0 
	P0_1480 = 0  
end 
λ_0p_980 = 980e-9 		# Comprimento de onda central da bomba em nm
λ_0p_1480 = 1480e-9 	# Comprimento de onda central da bomba em nm
Dλ_980 =1e-9 			# Largura de linha da bomba em nm ( 980 nm)
Dλ_1480 =1e-9 			# Largura de linha da bomba em nm (1480 nm)
λ_9801 =960e-9 			# nanometros - valor inferior λ para bombeamemnto
λ_9802 = 1000e-9 		# nanometros - valor superior λ para bombeamemnto 

""" end include("nicolau_Dados-M5-980.jl") """

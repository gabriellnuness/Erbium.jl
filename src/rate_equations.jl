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
function dPdz(P,β12, β21, n2, gama, Rho, λ)
  _ = P*(-β12*(1-n2) + β21*n2 - gama) + Rho/λ*β21*n2
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
Single pass configuration


Integrate the spectrum through erbium-doped optical fiber without mirrors in the setup,
unless the reflections from initial conditions.

 -----
|pump1|---------[wdm1]-----[erbium fiber]-----[wdm2]--------[OSA]
 -----		   /									\
              /										 |		-----
			 |										  -----|pump2|
			  ------[OSA]									-----

"""
function integrator_single_pass(input_parameters, fiber_parameters)
	
	h 	= 6.6e-34   # Planck's constant
	c 	= 3e8      	# Speed of light in Vacuum
	kB	= 1.38e-23  # Boltzmann constant

	P0_980_B 	= input_parameters.P0_980_B
	P0_980_F 	= input_parameters.P0_980_F
	P0_1480_B 	= input_parameters.P0_1480_B
	P0_1480_F 	= input_parameters.P0_1480_F

	RL_980 	= fiber_parameters.RL_980
	RL_1480 = fiber_parameters.RL_1480
	R0_980 	= fiber_parameters.R0_980
	R0_1480 = fiber_parameters.R0_1480
	T 		= fiber_parameters.T
	Ksigma 	= fiber_parameters.Ksigma
	Z		= fiber_parameters.Z
	ϵ1		= fiber_parameters.ϵ1
	ϵ2		= fiber_parameters.ϵ2
	gama 	= fiber_parameters.gama
	G		= fiber_parameters.G
	L		= fiber_parameters.L

	Rho0 	= Z*ϵ1*ϵ2


	z_intervals = 300
	λ_intervals = 300
	maximum_iterations = 500
	Nc = Int(λ_intervals/4)  # spectrum cut to ignore pump spectrum
	tol = 0.001

	λ1 = 1450e-9  # [m]
	λ2 = 1600e-9  # [m]

	λ0_pump_980 	= 980e-9 # [m]
	λ0_pump_1480 	= 1480e-9 # [m]
	Δλ_980 			= 1e-9 # [m]
	Δλ_1480 		= 1e-9 # [m]
	λ_9801 			= 960e-9 # [m]
	λ_9802 			= 1000e-9 # [m]
	
	λ_arr = zeros(λ_intervals+1)
	λ_980_arr = zeros(λ_intervals+1)
	β12 = zeros(λ_intervals+1)
	β21 = zeros(λ_intervals+1)
	β13 = zeros(λ_intervals+1)
	n2 = ones(z_intervals+1)
	z_arr = zeros(z_intervals+1)
	P1480_forward = zeros(λ_intervals+1, z_intervals+1)
	P1480_backward = zeros(λ_intervals+1, z_intervals+1)
	P980_forward = zeros(λ_intervals+1, z_intervals+1)
	P980_backward = zeros(λ_intervals+1, z_intervals+1)
	P_pump_980_B = zeros(λ_intervals+1)
	P_pump_980_F = zeros(λ_intervals+1)
	P_pump_1480_B = zeros(λ_intervals+1)
	P_pump_1480_F = zeros(λ_intervals+1)
	conv = ones(maximum_iterations+1)
	PRF = zeros(λ_intervals+1, z_intervals+1)   # auxiliar vector to use in convergence
	PRB = zeros(λ_intervals+1, z_intervals+1)   # auxiliar vector to use in convergence
	Ppz = zeros(λ_intervals+1)          		# pump power propagating on fiber
	gain = zeros(λ_intervals+1, z_intervals+1)

	dz = L/z_intervals 
	dλ = (λ2 - λ1)/λ_intervals
	dλ_980= (λ_9802 - λ_9801)/λ_intervals

	# read absorption spectra 
	# TODO: implement other fibers using the optical_fiber.jl script
	package_path = dirname(dirname(pathof(Erbium)))
	data_fiber = readdlm("$(package_path)/data/M5_abs.txt", ',') 
	β13_Interpolate=linear_interpolation(data_fiber[1:580,1], data_fiber[1:580,2])
	β12_Interpolate = linear_interpolation(data_fiber[581:end,1], data_fiber[581:end,2])

	for λ=1:λ_intervals+1
		λ_arr[λ] = λ1 + (λ-1)*dλ
		λ_980_arr[λ] = λ_9801+(λ-1)*dλ_980
		β12[λ] = G*0.2303*β12_Interpolate(λ_arr[λ]*1e9)
		β13[λ] = G*0.2303*β13_Interpolate(λ_980_arr[λ]*1e9)
		
		P_pump_980_F[λ] = 2/Δλ_980*sqrt(log(2)/π)*P0_980_F*exp(-(4*log(2)*((λ_980_arr[λ]-λ0_pump_980)/Δλ_980)^2))
		P_pump_980_B[λ] = 2/Δλ_980*sqrt(log(2)/π)*P0_980_B*exp(-(4*log(2)*((λ_980_arr[λ]-λ0_pump_980)/Δλ_980)^2))
		P_pump_1480_F[λ] = 2/Δλ_1480*sqrt(log(2)/π)*P0_1480_F*exp(-(4*log(2)*((λ_arr[λ]-λ0_pump_1480)/Δλ_1480)^2))
		P_pump_1480_B[λ] = 2/Δλ_1480*sqrt(log(2)/π)*P0_1480_B*exp(-(4*log(2)*((λ_arr[λ]-λ0_pump_1480)/Δλ_1480)^2))
	end

	# Calculate β21 using McCumber's approximation
	for λ=1:λ_intervals+1
		β21[λ]=β12[λ]*exp(-(h*c)/(λ_arr[λ]*kB*T))
	end

	β12_Max=maximum(β12)
	β21_Max=maximum(β21)

	for λ=1:λ_intervals+1
		β21[λ]=Ksigma*β12_Max/β21_Max*β21[λ]
	end

	Rho=Rho0/trapz(β21,dλ)

	# initialize variables
	for λ=1:λ_intervals+1     
		for z=1:z_intervals+1
			z_arr[z]=(z-1)*dz
			P1480_forward[λ,z] = 0
			P1480_backward[λ,z]= 0
			P980_forward[λ,z]=0;
			P980_backward[λ,z]=0;
			n2[z] = 0    
		end
	end

	global Refa=1
	global Refb=1
	conv[1]=1

	# begin integration of single_pass()
	for k = 1:maximum_iterations
		Refa, Refb
		# Propagation ---------------------------------------------->>>>
	
		# contourn conditions
		for λ=1:λ_intervals+1
			P1480_forward[λ,1]=R0_1480*P1480_backward[λ,1]+P_pump_1480_F[λ]
			P980_forward[λ,1]=R0_980*P980_backward[λ,1]+P_pump_980_F[λ]   
		end # λ

		for z=1:z_intervals+1
			global w12=0
			global w21=0
			global w13=0
			for λ=1:λ_intervals+1
				w12=w12+(P1480_forward[λ,z]+P1480_backward[λ,z])*λ_arr[λ]*β12[λ]*dλ
				w21=w21+(P1480_forward[λ,z]+P1480_backward[λ,z])*λ_arr[λ]*β21[λ]*dλ
				w13=w13+(P980_forward[λ,z]+P980_backward[λ,z])*λ_980_arr[λ]*β13[λ]*dλ_980
			end #λ
			n2[z]=(w13+w12)/((w13+w12+w21)+Z)
		end #z

		for z=1:z_intervals
			for λ=1:λ_intervals+1
				fp1 = dPdz(P980_forward[λ,z],		  β13[λ], 0, n2[z], gama, 0, λ_980_arr[λ])*dz
				fp2 = dPdz(P980_forward[λ,z]+fp1/2, β13[λ], 0, n2[z], gama, 0, λ_980_arr[λ])*dz
				fp3 = dPdz(P980_forward[λ,z]+fp2/2, β13[λ], 0, n2[z], gama, 0, λ_980_arr[λ])*dz
				fp4 = dPdz(P980_forward[λ,z]+fp3,	  β13[λ], 0, n2[z], gama, 0, λ_980_arr[λ])*dz
				P980_forward[λ,z+1] = P980_forward[λ,z] + (fp1+ 2*fp2 + 2*fp3 + fp4)/6
			end
			Ppz[z] = maximum(P980_forward[:,z])
			Ppz[z_intervals+1] = maximum(P980_forward[:,z_intervals+1])

	end

		for z=1:z_intervals          
			for λ=1:λ_intervals+1
				fp1 = dPdz(P1480_forward[λ,z],β12[λ],β21[λ],n2[z],gama,Rho, λ_arr[λ])*dz
				fp2 = dPdz(P1480_forward[λ,z]+fp1/2,β12[λ],β21[λ],n2[z],gama,Rho, λ_arr[λ])*dz
				fp3 = dPdz(P1480_forward[λ,z]+fp2/2,β12[λ],β21[λ],n2[z],gama,Rho,λ_arr[λ])*dz
				fp4 = dPdz(P1480_forward[λ,z]+fp3, β12[λ],β21[λ],n2[z],gama,Rho,λ_arr[λ])*dz
				P1480_forward[λ,z+1] = P1480_forward[λ,z]+(fp1+ 2*fp2 + 2*fp3 + fp4)/6
			end
		end 

		
	
		# <<<<---------------------------------------------- Propagation 

		for λ=1:λ_intervals+1
		# contourn conditions
			P1480_backward[λ,z_intervals+1]=RL_1480*P1480_forward[λ,z_intervals+1]+P_pump_1480_B[λ]
			P980_backward[λ,z_intervals+1]=RL_980*P980_forward[λ,z_intervals+1]+P_pump_980_B[λ]
		end # λ

		for z=1:z_intervals+1
			global w12=0
			global w21=0
			global w13=0
			for λ=1:λ_intervals+1
				w12=w12+(P1480_forward[λ,z]+P1480_backward[λ,z])*λ_arr[λ]*β12[λ]*dλ
				w21=w21+(P1480_forward[λ,z]+P1480_backward[λ,z])*λ_arr[λ]*β21[λ]*dλ
				w13=w13+(P980_forward[λ,z]+P980_backward[λ,z])*λ_980_arr[λ]*β13[λ]*dλ_980
			end #λ
			n2[z]=(w13+w12)/((w13+w12+w21)+Z)
		end #z

		for z=1:z_intervals
			jj=z_intervals+1-z          
			for λ=1:λ_intervals+1
				fp1 = dPdz(P980_backward[λ,jj],β13[λ], 0,n2[jj],gama,0, λ_980_arr[λ])*dz
				fp2 = dPdz(P980_backward[λ,jj]+fp1/2,β13[λ], 0,n2[jj],gama,0, λ_980_arr[λ])*dz
				fp3 = dPdz(P980_backward[λ,jj]+fp2/2,β13[λ], 0,n2[jj],gama,0,λ_980_arr[λ])*dz
				fp4 = dPdz(P980_backward[λ,jj]+fp3, β13[λ], 0,n2[jj],gama,0,λ_980_arr[λ])*dz
				P980_backward[λ,jj]=P980_backward[λ,jj+1]+(fp1+ 2*fp2 + 2*fp3 + fp4)/6
			end
		end 


		for z=1:z_intervals    
			jj=z_intervals+1-z
			for λ=1:λ_intervals+1
				fp1 = dPdz(P1480_backward[λ,jj+1],β12[λ],β21[λ],n2[jj+1],gama,Rho,λ_arr[λ])*dz
				fp2 = dPdz(P1480_backward[λ,jj+1]+fp1/2,β12[λ],β21[λ],n2[jj+1],gama,Rho,λ_arr[λ])*dz
				fp3 = dPdz(P1480_backward[λ,jj+1]+fp2/2,β12[λ],β21[λ],n2[jj+1],gama,Rho,λ_arr[λ])*dz
				fp4 = dPdz(P1480_backward[λ,jj+1]+fp3,β12[λ],β21[λ],n2[jj+1],gama,Rho,λ_arr[λ])*dz
				P1480_backward[λ,jj]=P1480_backward[λ,jj+1]+(fp1 + 2*fp2 + 2*fp3 + fp4)/6
			end
		end #z     

		global Refb=(maximum(P1480_forward[Nc:λ_intervals+1,z_intervals+1]))
		Max=Refb
		conv[k+1]=abs((Refb-Refa)/Refa)
		if conv[k+1] > 2
			conv[k+1]=2
		end

		println("k = " ,k,"  conv = ",conv[k+1])
		if k >6
			if  (conv[k-1]< tol)&&(conv[k]< tol)&&(conv[k+1]<tol)
				break
			end
		end     
		global Refa=Refb

		if k>20
			for λ = 1:λ_intervals+1
				for z=1:z_intervals+1
					P1480_forward[λ,z]=(P1480_forward[λ,z]+PRF[λ,z])/2
					P1480_backward[λ,z]=(P1480_backward[λ,z]+PRB[λ,z])/2
				end
			end
		end     

		for λ = 1:λ_intervals+1
			for z=1:z_intervals+1
				PRF[λ,z]=P1480_forward[λ,z]
				PRB[λ,z]=P1480_backward[λ,z]
			end
		end

	end 

	for λ=1:λ_intervals+1
		for z=1:z_intervals+1
			gain[λ,z]=L*(β21[λ]*n2[z]-β12[λ]*(1-n2[z]))
		end
	end


	return (power_forward = P1480_forward[:,z_intervals+1], 		# return named tuple
			power_backward = P1480_backward[:,1],
			λ = λ_arr,
			gain = gain)
end # end test_integrator_single_pass()



function calc_spec_parameters(power_forward, power_backward, λ)

	c 	= 3e8      	# speed of light in vacuum

	dλ = λ[2]-λ[1]
	ITf = trapz(power_forward,dλ)  					# ∫ P1480_forward(λ) dλ 
	ITb = trapz(power_backward,dλ)     				# ∫ P1480_backward(λ) dλ 
	ITlf = trapz(power_forward.*λ,dλ) 				# ∫ λ P1480_forward(λ) dλ 
	ITlb = trapz(power_backward.*λ,dλ) 				# ∫ λ P1480_backward(λ) dλ 
	ITf2 = trapz(power_forward.*power_forward,dλ) 	# ∫ [P1480_forward(λ)]² dλ 
	ITb2 = trapz(power_backward.*power_backward,dλ) # ∫ [P1480_backward(λ)]² dλ 

	Potf = ITf  				# Potência saída Forward
	Potb = ITb   				# Potência saída Backward
	λmf = ITlf/ITf  			# <λ> Forward
	λmb = ITlb/ITb  			# <λ> Backward
	DλEff_f = ITf^2/ITf2 		# Δλeff Forward
	DλEff_b = ITb^2/ITb2 		# Δλeff Backward
	RIN_f = λmf^2 /(c*DλEff_f)
	RIN_b = λmb^2 /(c*DλEff_b)
	
	return (Potf=Potf, Potb=Potb, λmf=λmf, λmb=λmb, DλEff_f=DλEff_f, DλEff_b=DλEff_b, RIN_f=RIN_f, RIN_b=RIN_b)
end


function trapz(A,dl)
    S=0
    for i=1:size(A,1)-1
         S=S+(A[i]+A[i+1])/2*dl
    end
    return S   
end

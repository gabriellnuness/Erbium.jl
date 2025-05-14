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
function integrate_single_pass()

	###Begin: Definição de variáveis
	λ=zeros(Nl+1)       # λ na região de 1480 nm
	λ_980=zeros(Nl+1)   # λ na região de 980 nm
	β12=zeros(Nl+1)     # Coeficiente de absorção na região de 1480 nm
	β21=zeros(Nl+1)     # Coeficiente de emissão na região de 1480 nm
	β13=zeros(Nl+1)   # Coeficiente de absorção na região de 980 nm
	wdm=zeros(Nl+1)     # Função transmissão do WDM
	n2=ones(Nz+1)       # População normalizada do nível 2
	z=zeros(Nz+1)       # Posição ao longo da propagação
	P1480F=zeros(Nl+1, Nz+1)     # Potência que se propaga para a direita
	P1480B=zeros(Nl+1, Nz+1)     # Potência que se propaga para a esquerda
	P980F=zeros(Nl+1, Nz+1)    # Potência de bomba que se propaga para a direita
	P980B=zeros(Nl+1, Nz+1)    # Potência de bomba que se propaga para a esquerda
	P_980=zeros(Nl+1)        # Potência de bomba total em 980 nm
	P_1480=zeros(Nl+1)       # Potência de bomba total em 1480 nm
	conv=ones(Nk+1)          # Verificador de convergência
	PRF=zeros(Nl+1, Nz+1)    # Vetor auxiliar para acelerar convergência
	PRB=zeros(Nl+1, Nz+1)    # Vetor auxiliar para acelerar convergência
	Ppz=zeros(Nl+1)          # Potência de bomba que se propaga na fibra
	gain=zeros(Nl+1, Nz+1)   # Distribuição de ganho (λ,z)
	###End: Definição de variáveis

	dz = L/Nz # Número de incrementos em comprimento; 
	dλ = (λ2 - λ1)/Nl #nanometros - Incremento em λ
	dλ_980= (λ_9802 - λ_9801)/Nl #nanometros - Incremento em λ - bomba

	###Begin: Leituras dos espectros de Absorção, Emissão estimulada e de Bomba
	data_Ab = readdlm("dual_pump/M-5 abs.txt") 
	β12_Interpolate=linear_interpolation(data_Ab[:,1], data_Ab[:,2])
	
	data_P = readdlm("dual_pump/M-5 980.txt")
	β13_Interpolate = linear_interpolation(data_P[:,1], data_P[:,2])


	# Os contadores serão associados a 
		# i ==> λ
		# j ==> z
		# k ==> Perações

	#####Begin: Definir densidades espectrais #####
	for i=1:Nl+1
		λ[i]=λ1 + (i-1)*dλ #Discretização do comprimento de onda
		λ_980[i]=λ_9801+(i-1)*dλ_980 #Discretização do comprimento de onda de bomba
		β12[i]=G*0.2303*β12_Interpolate(λ[i]*1e9) #βa
		β13[i]=G*0.2303*β13_Interpolate(λ_980[i]*1e9) #βa de bomba
		
		### Definição da distribuição espectral de Pensidade da bomba 
		### na entrada da fibra
		P_980[i] =2/Dλ_980*sqrt(log(2)/π)*P0_980*exp(-(4*log(2)*((λ_980[i]-λ_0p_980)/Dλ_980)^2))
		P_1480[i] =2/Dλ_1480*sqrt(log(2)/π)*P0_1480*exp(-(4*log(2)*((λ[i]-λ_0p_1480)/Dλ_1480)^2))
	end
	#####End: Definir densidades espectrais #####

	#####Begin: Calcular β21 usando aproximação de McCumber ########
		for i=1:Nl+1
			β21[i]=β12[i]*exp(-(h*c)/(λ[i]*kB*T))
		end

		β12_Max=maximum(β12)
		β21_Max=maximum(β21)

		for i=1:Nl+1
			β21[i]=Ksigma*β12_Max/β21_Max*β21[i]
		end
	#####End: Calcular β21 usando aproximação de McCumber ########
	SS=Integral(β21,dλ)
	Rho=Rho0/SS

	###Begin: Atribuição de valores iniciais ####
	for i=1:Nl+1     
		for j=1:Nz+1
			z[j]=(j-1)*dz  #Atribuição dos valores de z
			P1480F[i,j] = 0 #Valores iniciais da Pensidade forw12rd
			P1480B[i,j]= 0 #Valores iniciais da Pensidade backw12rd
			P980F[i,j]=0;
			P980B[i,j]=0;
			n2[j] = 0  # Atribuicao valores iniciais de N2     
		end
	end  
	###End: Atribuição de valores iniciais ####

	##########Begin: Resolução da equação diferencial  ############

	global Refa=1
	global Refb=1
	conv[1]=1

	for k = 1:Nk
		Refa, Refb
		# Propagação para a direita
	
		#Condições de contorno
		for i=1:Nl+1
			P1480F[i,1]=R0_1480*P1480B[i,1]#+P_1480[i]
			P980F[i,1]=R0_980*P980B[i,1]+P_980[i]   
		end # i

		for j=1:Nz+1
			global w12=0
			global w21=0
			global w13=0
			for i=1:Nl+1
				w12=w12+(P1480F[i,j]+P1480B[i,j])*λ[i]*β12[i]*dλ
				w21=w21+(P1480F[i,j]+P1480B[i,j])*λ[i]*β21[i]*dλ
				w13=w13+(P980F[i,j]+P980B[i,j])*λ_980[i]*β13[i]*dλ_980
			end #i
			n2[j]=(w13+w12)/((w13+w12+w21)+Z)
		end #j

		for j=1:Nz
			for i=1:Nl+1
				fp1 = dPdz(P980F[i,j],β13[i],0,n2[j],gama,0, λ_980[i])*dz
				fp2 = dPdz(P980F[i,j]+fp1/2,β13[i], 0,n2[j],gama,0, λ_980[i])*dz
				fp3 = dPdz(P980F[i,j]+fp2/2,β13[i], 0,n2[j],gama,0,λ_980[i])*dz
				fp4 = dPdz(P980F[i,j]+fp3, β13[i],0,n2[j],gama,0,λ_980[i])*dz
				P980F[i,j+1]=P980F[i,j]+(fp1+ 2*fp2 + 2*fp3 + fp4)/6
			end
			Ppz[j] = maximum(P980F[:,j])
			Ppz[Nz+1] = maximum(P980F[:,Nz+1])

	end

		for j=1:Nz          
			for i=1:Nl+1
				fp1 = dPdz(P1480F[i,j],β12[i],β21[i],n2[j],gama,Rho, λ[i])*dz
				fp2 = dPdz(P1480F[i,j]+fp1/2,β12[i],β21[i],n2[j],gama,Rho, λ[i])*dz
				fp3 = dPdz(P1480F[i,j]+fp2/2,β12[i],β21[i],n2[j],gama,Rho,λ[i])*dz
				fp4 = dPdz(P1480F[i,j]+fp3, β12[i],β21[i],n2[j],gama,Rho,λ[i])*dz
				P1480F[i,j+1] = P1480F[i,j]+(fp1+ 2*fp2 + 2*fp3 + fp4)/6
			end
		end 

		
	#Propagação para a esquerda
	
		#Condições de contorno
		for i=1:Nl+1
			P1480B[i,Nz+1]=RL_1480*P1480F[i,Nz+1]+P_1480[i]
			P980B[i,Nz+1]=RL_980*P980F[i,Nz+1]#+P_980[i]
		end # i

		for j=1:Nz+1
			global w12=0
			global w21=0
			global w13=0
			for i=1:Nl+1
				w12=w12+(P1480F[i,j]+P1480B[i,j])*λ[i]*β12[i]*dλ
				w21=w21+(P1480F[i,j]+P1480B[i,j])*λ[i]*β21[i]*dλ
				w13=w13+(P980F[i,j]+P980B[i,j])*λ_980[i]*β13[i]*dλ_980
			end #i
			n2[j]=(w13+w12)/((w13+w12+w21)+Z)
		end #j

	for j=1:Nz
		jj=Nz+1-j          
		for i=1:Nl+1
			fp1 = dPdz(P980B[i,jj],β13[i], 0,n2[jj],gama,0, λ_980[i])*dz
			fp2 = dPdz(P980B[i,jj]+fp1/2,β13[i], 0,n2[jj],gama,0, λ_980[i])*dz
			fp3 = dPdz(P980B[i,jj]+fp2/2,β13[i], 0,n2[jj],gama,0,λ_980[i])*dz
			fp4 = dPdz(P980B[i,jj]+fp3, β13[i], 0,n2[jj],gama,0,λ_980[i])*dz
			P980B[i,jj]=P980B[i,jj+1]+(fp1+ 2*fp2 + 2*fp3 + fp4)/6
		end
	end 


	for j=1:Nz    
			jj=Nz+1-j
			for i=1:Nl+1
				fp1 = dPdz(P1480B[i,jj+1],β12[i],β21[i],n2[jj+1],gama,Rho,λ[i])*dz
				fp2 = dPdz(P1480B[i,jj+1]+fp1/2,β12[i],β21[i],n2[jj+1],gama,Rho,λ[i])*dz
				fp3 = dPdz(P1480B[i,jj+1]+fp2/2,β12[i],β21[i],n2[jj+1],gama,Rho,λ[i])*dz
				fp4 = dPdz(P1480B[i,jj+1]+fp3,β12[i],β21[i],n2[jj+1],gama,Rho,λ[i])*dz
				P1480B[i,jj]=P1480B[i,jj+1]+(fp1 + 2*fp2 + 2*fp3 + fp4)/6
			end
		end #j     

		global Refb=(maximum(P1480F[Nc:Nl+1,Nz+1]))
		Max=Refb
		conv[k+1]=abs((Refb-Refa)/Refa)
		if conv[k+1] > 2
			conv[k+1]=2
		end

		println("k = " ,k,"  conv = ",conv[k+1])
		if k >6
			if  (conv[k-1]< Tol)&&(conv[k]< Tol)&&(conv[k+1]< Tol)
				break
			end
		end     
	global Refa=Refb

	if k>20
		for i = 1:Nl+1
			for j=1:Nz+1
				P1480F[i,j]=(P1480F[i,j]+PRF[i,j])/2
				P1480B[i,j]=(P1480B[i,j]+PRB[i,j])/2
			end
		end
	end     


	#plot(P1480F[:,Nz+1])
	#plot(P1480B[:,1])

		for i = 1:Nl+1
			for j=1:Nz+1
				PRF[i,j]=P1480F[i,j]
				PRB[i,j]=P1480B[i,j]
			end
		end

	end #Fim loop k  


end
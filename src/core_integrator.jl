#= 
Simulação de uma fonte de ASE usando fibra dopada com érbio:
- Bombeamento com λ = 980 nm e/ou λ = 1480 nm;
- Cálculo do espectro de emissão dependente da Temperatura
  usando aproximação de McCumber;
=#


#include("Dados-M12.jl")

function dPdz(P,β12, β21, n2, gama, Rho, λ)
  _=P*(-β12*(1-n2)+β21*n2-gama)+Rho/λ*β21*n2
end

function Integral(A,dl)
    S=0
    for i=1:size(A,1)-1
         S=S+(A[i]+A[i+1])/2*dl
    end
    return S   
end





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
     data_fiber = readdlm("data/M5_abs.txt", ',') 
     β13_Interpolate=linear_interpolation(data_fiber[1:580,1], data_fiber[1:580,2])
     β12_Interpolate = linear_interpolation(data_fiber[581:end,1], data_fiber[581:end,2])


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
             fp1 = dPdz(P980F[i,j],		  β13[i], 0, n2[j], gama, 0, λ_980[i])*dz
             fp2 = dPdz(P980F[i,j]+fp1/2, β13[i], 0, n2[j], gama, 0, λ_980[i])*dz
             fp3 = dPdz(P980F[i,j]+fp2/2, β13[i], 0, n2[j], gama, 0, λ_980[i])*dz
             fp4 = dPdz(P980F[i,j]+fp3,	  β13[i], 0, n2[j], gama, 0, λ_980[i])*dz
             P980F[i,j+1] = P980F[i,j] + (fp1+ 2*fp2 + 2*fp3 + fp4)/6
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
	for i=1:Nl+1
	#Condições de contorno
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

end # end integrate_single_pass()

######### Até aqui foram calculadas P1480 e P980 Forward e Backward #######  

#####Begin: Cálculo da distribuição de ganho #####
for i=1:Nl+1
     for j=1:Nz+1
          gain[i,j]=L*(β21[i]*n2[j]-β12[i]*(1-n2[j]))
     end
end
#####End: Cálculo da distribuição de ganho #####

####Begin: Cálculo de Pot Saída, <λ> e Δλeff 
ITf=Integral(P1480F[Nc:end,Nz+1],dλ)  ### ∫ P1480F(λ) dλ ####
ITb=Integral(P1480B[Nc:end,1],dλ)     ### ∫ P1480B(λ) dλ ####
ITlf=Integral(P1480F[Nc:end,Nz+1].*λ[Nc:end],dλ) ### ∫ λ P1480F(λ) dλ ####
ITlb=Integral(P1480B[Nc:end,1].*λ[Nc:end],dλ) ### ∫ λ P1480B(λ) dλ ####
ITf2=Integral(P1480F[Nc:end,Nz+1].*P1480F[Nc:end,Nz+1],dλ) ### ∫ [P1480F(λ)]² dλ ####
ITb2=Integral(P1480B[Nc:end,1].*P1480B[Nc:end,1],dλ)  ### ∫ [P1480B(λ)]² dλ ####

Potf =ITf  ### Potência saída Forw12rd
Potb=ITb   ### Potência saída Backw12rd
λmf=ITlf/ITf  ### <λ> Forw12rd
λmb=ITlb/ITb  ### <λ> Backw12rd
DλEff_f= ITf^2/ITf2 ### Δλeff Forw12rd
DλEff_b= ITb^2/ITb2 ### Δλeff Backw12rd
RIN_f = λmf^2 /(c*DλEff_f)
RIN_b = λmb^2 /(c*DλEff_b)
####End: Cálculo de Pot Saída, <λ> e Δλeff

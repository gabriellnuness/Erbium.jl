

"""
	differential equation for 980 pumping scheme
"""
function dPdz(P,β_a, β_e, n2, gama, Rho, λ)
  _=P*(-β_a*(1-n2)+β_e*n2-gama)+Rho/λ*β_e*n2
end



"""
	trapezoidal integration?
"""
function Integral(A,dl)
    S=0
    for i = 1 : size(A,1)-1
         S = S + (A[i] + A[i+1]) / 2 * dl
    end
    return S   
end





###Begin: Definição de variáveis
     λ=zeros(Nλ+1)       # λ na região de 1480 nm
     λ_980=zeros(Nλ+1)   # λ na região de 980 nm
     β_a=zeros(Nλ+1)     # Coeficiente de absorção na região de 1480 nm
     β_e=zeros(Nλ+1)     # Coeficiente de emissão na região de 1480 nm
     β_980=zeros(Nλ+1)   # Coeficiente de absorção na região de 980 nm
     wdm=zeros(Nλ+1)     # Função transmissão do WDM
     n2=ones(Nz+1)       # População normalizada do nível 2
     z=zeros(Nz+1)       # Posição ao longo da propagação
     Pf=zeros(Nλ+1, Nz+1)     # Potência que se propaga para a direita
     Pb=zeros(Nλ+1, Nz+1)     # Potência que se propaga para a esquerda
     Ppf=zeros(Nλ+1, Nz+1)    # Potência de bomba que se propaga para a direita
     Ppb=zeros(Nλ+1, Nz+1)    # Potência de bomba que se propaga para a esquerda
     P_980=zeros(Nλ+1)        # Potência de bomba total em 980 nm
     P_1480=zeros(Nλ+1)       # Potência de bomba total em 1480 nm
     conv=ones(max_iter+1)          # Verificador de convergência
     PExp=zeros(Nλ+1)         # Potência experimental acoplada
     PRF=zeros(Nλ+1, Nz+1)    # Vetor auxiliar para acelerar convergência
     PRB=zeros(Nλ+1, Nz+1)    # Vetor auxiliar para acelerar convergência
     Ppz=zeros(Nλ+1)          # Potência de bomba que se propaga na fibra
     gain=zeros(Nλ+1, Nz+1)   # Distribuição de ganho (λ,z)
###End: Definição de variáveis

dz = L / Nz 					# Número de incrementos em comprimento; 
dλ = (λ2 - λ1)/Nλ 				# nanometros - Incremento em λ
dλ_980= (λ_9802 - λ_9801)/Nλ 	# nanometros - Incremento em λ - bomba



###Begin: Leituras dos espectros de Absorção, Emissão estimulada e de Bomba
    #  data_Ab = readdlm(Dir*"/M-5 abs.txt") 
    #  β_a_Interpolate=linear_interpolation(data_Ab[:,1], data_Ab[:,2])
     
    #  data_P = readdlm(Dir*"/M-5 980.txt")
    #  β_980_Interpolate = linear_interpolation(data_P[:,1], data_P[:,2])


fiber_model = "m5"
(β_abs0,β_emis0,τ21,τ3,fiber_diameter,NA,total_population,η) = optical_fiber(fiber_model)

figure(1)
    plot(β_abs0[:,1], β_abs0[:,2], label="$fiber_model absorption")
    plot(β_emis0[:,1], β_emis0[:,2], label="$fiber_model emission")
	ylabel(L"$\beta$ [dB/m]")
	xlabel("wavelength [m]")
	legend()

β_a_Interpolate = linear_interpolation(β_abs0[:,1], β_abs0[:,2])
β_980_Interpolate = linear_interpolation(β_abs0[:,1], β_abs0[:,2])
# β_980_Interpolate = linear_interpolation(β_emis0[:,1], β_emis0[:,2])



####Begin: Leitura resultado Exp ######

arq_Exp_F = readdlm(FileExp_F,',')
arq_Exp_B = readdlm(FileExp_B,',')
data_F = arq_Exp_F[2:end, 2]
data_B = arq_Exp_B[2:end, 2]
λ_Exp_F = arq_Exp_F[2:end, 1]
#λ_Exp_B=arq_Exp_B[:,1]
# dλ_Exp_F= λ_Exp_F[2]-λ_Exp_F[1]
#dλ_Exp_B= λ_Exp_B[2]-λ_Exp_B[1]

####End: Leitura resultado Exp ######

figure(2)
	plot(λ_Exp_F, data_F, label="forward")
	plot(λ_Exp_F, data_B, label="backward")
	legend()
	ylabel("experimental spectrum [dBm]")
	xlabel("wavelength [m]")

# Os contadores serão associados a 
     # i ==> λ
     # j ==> z
     # k ==> Perações

#####Begin: Definir densidades espectrais #####
for i = 1 : Nλ+1
     λ[i] = λ1 + (i-1)*dλ 							# Discretização do comprimento de onda
     λ_980[i] = λ_9801+(i-1)*dλ_980 				# Discretização do comprimento de onda de bomba
    #  β_a[i] = G*0.2303*β_a_Interpolate(λ[i]*1e9)	# β absorption
     β_a[i] = 0.2303*β_a_Interpolate(λ[i])	# β absorption
    #  β_980[i] = G*0.2303*β_980_Interpolate(λ_980[i]*1e9) #βabsorption pump laser
     β_980[i] = G*0.2303*β_980_Interpolate(λ_980[i]) #βabsorption pump laser
     
     ### Definição da distribuição espectral de Pensidade da bomba 
     ### na entrada da fibra
     P_980[i] = 2/Dλ_980*sqrt(log(2)/π)*P0_980*exp(-(4*log(2)*((λ_980[i]-λ_0p_980)/Dλ_980)^2))
     P_1480[i] = 2/Dλ_1480*sqrt(log(2)/π)*P0_1480*exp(-(4*log(2)*((λ[i]-λ_0p_1480)/Dλ_1480)^2))
end
#####End: Definir densidades espectrais #####

figure()
	plot(λ, β_a)
	plot(λ_980, β_980, color=custom_plot_colors[1], label="M5 absorption separated")

figure()
	plot(λ_980, P_980, label="Power of 980 nm")
	# plot(λ, P_1480, label="Power of 1480 nm")
	legend()





"""
	Begin: Calcular β_e usando aproximação de McCumber
"""

for i=1:Nλ+1
	β_e[i] = β_a[i]*exp(-(h*c)/(λ[i]*kB*T))
end

β_a_Max=maximum(β_a)
β_e_Max=maximum(β_e)

for i=1:Nλ+1
	β_e[i]=Ksigma*β_a_Max/β_e_Max*β_e[i]
end

figure()
	plot(λ, β_a, label="given M5 absorption")
	plot(λ, β_e, label="calculates M5 emission")
	legend()
	ylabel("dB/m")
	xlabel("wavelength [m]")

figure()
	plot(λ, β_a, label="given M5 absorption")
	plot(λ, β_e, label="calculates M5 emission")
 	plot(β_abs0[:,1], 0.2303*β_abs0[:,2], label="$fiber_model absorption", color=custom_plot_colors[1], "--")
    plot(β_emis0[:,1], 0.2303*β_emis0[:,2], label="$fiber_model emission", color=custom_plot_colors[2], "--")
	legend()
	ylabel("dB/m")
	xlabel("wavelength [m]")

"""
	end: Calcular β_e usando aproximação de McCumber
"""


SS = Integral(β_e, dλ)
Rho = Rho0 / SS

###Begin: Atribuição de valores iniciais ####
for i=1:Nλ+1     
     for j = 1:Nz+1
          z[j]=(j-1)*dz  #Atribuição dos valores de z
          Pf[i,j] = 0 #Valores iniciais da Pensidade forward
          Pb[i,j]= 0 #Valores iniciais da Pensidade backward
          Ppf[i,j]=0;
          Ppb[i,j]=0;
          n2[j] = 0  # Atribuicao valores iniciais de N2     
     end
end  
###End: Atribuição de valores iniciais ####

##########Begin: Resolução da equação diferencial  ############

global Refa=1
global Refb=1
conv[1]=1

for k = 1:max_iter
     Refa, Refb
     # Propagação para a direita
  
     #Condições iniciais
     for i=1:Nλ+1
          Pf[i,1]=R0*Pb[i,1]+P_1480[i]
          Ppf[i,1]=R0*Ppb[i,1]+P_980[i]   # <==== Ten umn erro aqui, λs diferentes, 
                                             #mas que não deve afetar o caso de λ = 14980 n
     end # i

     for j=1:Nz+1
          global wa=0
          global we=0
          global wp=0
          for i=1:Nλ+1
               wa=wa+(Pf[i,j]+Pb[i,j])*λ[i]*β_a[i]*dλ
               we=we+(Pf[i,j]+Pb[i,j])*λ[i]*β_e[i]*dλ
               wp=wp+(Ppf[i,j]+Ppb[i,j])*λ_980[i]*β_980[i]*dλ_980
          end #i
          n2[j]=(wp+wa)/((wp+wa+we)+Z)
     end #j

     for j=1:Nz
        for i=1:Nλ+1
             fp1 = dPdz(Ppf[i,j],β_980[i],0,n2[j],gama,0, λ_980[i])*dz
             fp2 = dPdz(Ppf[i,j]+fp1/2,β_980[i], 0,n2[j],gama,0, λ_980[i])*dz
             fp3 = dPdz(Ppf[i,j]+fp2/2,β_980[i], 0,n2[j],gama,0,λ_980[i])*dz
             fp4 = dPdz(Ppf[i,j]+fp3, β_980[i],0,n2[j],gama,0,λ_980[i])*dz
             Ppf[i,j+1]=Ppf[i,j]+(fp1+ 2*fp2 + 2*fp3 + fp4)/6
        end
        Ppz[j] = maximum(Ppf[:,j+1])
   end

     for j=1:Nz          
          for i=1:Nλ+1
               fp1 = dPdz(Pf[i,j],β_a[i],β_e[i],n2[j],gama,Rho, λ[i])*dz
               fp2 = dPdz(Pf[i,j]+fp1/2,β_a[i],β_e[i],n2[j],gama,Rho, λ[i])*dz
               fp3 = dPdz(Pf[i,j]+fp2/2,β_a[i],β_e[i],n2[j],gama,Rho,λ[i])*dz
               fp4 = dPdz(Pf[i,j]+fp3, β_a[i],β_e[i],n2[j],gama,Rho,λ[i])*dz
               Pf[i,j+1] = Pf[i,j]+(fp1+ 2*fp2 + 2*fp3 + fp4)/6
          end
     end 

       
#Propagação para a esquerda
   
     #Condições iniciais
     for i=1:Nλ+1
          Pb[i,Nz+1]=RL*Pf[i,Nz+1]
          Ppb[i,Nz+1]=RL*Ppf[i,Nz+1]
     end # i

     for j=1:Nz+1
          global wa=0
          global we=0
          global wp=0
          for i=1:Nλ+1
               wa=wa+(Pf[i,j]+Pb[i,j])*λ[i]*β_a[i]*dλ
               we=we+(Pf[i,j]+Pb[i,j])*λ[i]*β_e[i]*dλ
               wp=wp+(Ppf[i,j]+Ppb[i,j])*λ_980[i]*β_980[i]*dλ_980
          end #i
          n2[j]=(wp+wa)/((wp+wa+we)+Z)
     end #j

for j=1:Nz
     jj=Nz+1-j          
    for i=1:Nλ+1
         fp1 = dPdz(Ppb[i,jj],β_980[i], 0,n2[jj],gama,0, λ_980[i])*dz
         fp2 = dPdz(Ppb[i,jj]+fp1/2,β_980[i], 0,n2[jj],gama,0, λ_980[i])*dz
         fp3 = dPdz(Ppb[i,jj]+fp2/2,β_980[i], 0,n2[jj],gama,0,λ_980[i])*dz
         fp4 = dPdz(Ppb[i,jj]+fp3, β_980[i], 0,n2[jj],gama,0,λ_980[i])*dz
         Ppb[i,jj]=Ppb[i,jj+1]+(fp1+ 2*fp2 + 2*fp3 + fp4)/6
    end
end 


for j=1:Nz    
          jj=Nz+1-j
          for i=1:Nλ+1
               fp1 = dPdz(Pb[i,jj+1],β_a[i],β_e[i],n2[jj+1],gama,Rho,λ[i])*dz
               fp2 = dPdz(Pb[i,jj+1]+fp1/2,β_a[i],β_e[i],n2[jj+1],gama,Rho,λ[i])*dz
               fp3 = dPdz(Pb[i,jj+1]+fp2/2,β_a[i],β_e[i],n2[jj+1],gama,Rho,λ[i])*dz
               fp4 = dPdz(Pb[i,jj+1]+fp3,β_a[i],β_e[i],n2[jj+1],gama,Rho,λ[i])*dz
               Pb[i,jj]=Pb[i,jj+1]+(fp1 + 2*fp2 + 2*fp3 + fp4)/6
          end
     end #j     

     global Refb=(maximum(Pf[Nc:Nλ+1,Nz+1]))
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
     for i = 1:Nλ+1
          for j=1:Nz+1
               Pf[i,j]=(Pf[i,j]+PRF[i,j])/2
               Pb[i,j]=(Pb[i,j]+PRB[i,j])/2
          end
     end
end     


#plot(Pf[:,Nz+1])
#plot(Pb[:,1])

     for i = 1:Nλ+1
          for j=1:Nz+1
               PRF[i,j]=Pf[i,j]
               PRB[i,j]=Pb[i,j]
          end
     end

end #Fim loop k  

#########  End:  Resolução da equação diferencial #########
######### Até aqui foram calculadas Pforward e Pbackward #######  

#####Begin: Cálculo da distribuição de ganho #####
for i=1:Nλ+1
     for j=1:Nz+1
          gain[i,j]=L*(β_e[i]*n2[j]-β_a[i]*(1-n2[j]))
     end
end
#####End: Cálculo da distribuição de ganho #####

####Begin: Cálculo de Pot Saída, <λ> e Δλeff 
ITf=Integral(Pf[Nc:end,Nz+1],dλ)  ### ∫ PF(λ) dλ ####
ITb=Integral(Pb[Nc:end,1],dλ)     ### ∫ PB(λ) dλ ####
ITlf=Integral(Pf[Nc:end,Nz+1].*λ[Nc:end],dλ) ### ∫ λ PF(λ) dλ ####
ITlb=Integral(Pb[Nc:end,1].*λ[Nc:end],dλ) ### ∫ λ PB(λ) dλ ####
ITf2=Integral(Pf[Nc:end,Nz+1].*Pf[Nc:end,Nz+1],dλ) ### ∫ [PF(λ)]² dλ ####
ITb2=Integral(Pb[Nc:end,1].*Pb[Nc:end,1],dλ)  ### ∫ [PB(λ)]² dλ ####

Potf =ITf  ### Potência saída Forward
Potb=ITb   ### Potência saída Backward
λmf=ITlf/ITf  ### <λ> Forward
λmb=ITlb/ITb  ### <λ> Backward
DλEff_f= ITf^2/ITf2 ### Δλeff Forward
DλEff_b= ITb^2/ITb2 ### Δλeff Backward
RIN_f = λmf^2 /(c*DλEff_f)
RIN_b = λmb^2 /(c*DλEff_b)
####End: Cálculo de Pot Saída, <λ> e Δλeff

######Begin: Calcular a potência acoplada em dBm e experimental em Watts ###########
NRF=size(data_F,1); 
NRB=size(data_B,1)
  ### Número de pontos experimentais
PotExpF = zeros(NRF)     #### Pot experimental em watts
PotExpB = zeros(NRF)     #### Pot experimental em watts
PotdBmf = zeros(Nλ+1)  #### Forward
PotdBmb = zeros(Nλ+1)  #### Backward
PotCalcF = zeros(Nλ+1);
PotCalcB = zeros(Nλ+1);

for i=1:Nλ+1
     PotCalcF[i]=Pf[i,Nz+1]*1e3;  # Converter em mW
     PotCalcB[i]=Pb[i,1]*1e3;     # Converter em mW
     PotdBmf[i]=10*log10(Pf[i,Nz+1]*1e3)
     PotdBmb[i]=10*log10(Pb[i,1]*1e3)
 end

 for i=1:NRF
     PotExpF[i]= 1e-3*10^(data_F[i]/10);
end

for i=1:NRB
     PotExpB[i]= 1e-3*10^(data_B[i]/10);
end

######End: Calcular a potência acoplada em dBm ###########


#Normalização para comparação
MaxPdBmf=maximum(PotdBmf[Nc:end])
MaxPdBmb=maximum(PotdBmb[Nc:end])
MaxPExp_F=maximum(PotExpF[Nexp:end]);
MaxPExp_B=maximum(PotExpB[Nexp:end]);
MaxPotCalcF=maximum(PotCalcF[Nc:end]);
MaxPotCalcB=maximum(PotCalcB[Nc:end]);
MR=maximum(data_F)
MI=maximum(PotdBmb)


figure()
	title("experiment")
	plot(λ_Exp_F*1e9, PotExpF/MaxPExp_F, label="forward power normalized")
	plot(λ_Exp_F*1e9, PotExpB/MaxPExp_B, label="backward power normalized")
	xlim(1500, 1580)
	yscale("linear")


figure()
	title("model")
	plot(λ*1e9, PotCalcF / MaxPotCalcF, label="forward power normalized")
	plot(λ*1e9, PotCalcB / MaxPotCalcB, label="backward power normalized")
	xlim(1500, 1580)


# ######Begin: Escrever os resultados calculados na tela e em arquivo ######
# F = DataFrame()

#      F.Parametros = ["I_980 [mA] = ";
#           "P0 980 nm [mW] =  ";
#      "λ_0p [nm]=  ";
#      "Dλ_980  [nm] =  ";
#      "I-1480 [mA]  = "; 
#       "P0 1480 nm [mW] =  ";
#      "λ_0p [nm]=  ";
#      "Dλ_1480  [nm] =  ";
 
#      "RL =  ";
#      "R0 =  ";
#      "K0 = "
#      "Rho0 = "
#       "T [K] =  ";
#       "L [m] =  ";
#      "Pot Forward  [mW]  = "
#      "Pot backward [mW]  =  ";
#      "λ médio forward [nm] =  ";
#      "λ médio backward [nm] = ";
#      "DλEff_f  [nm] = ";
#      "DλEff_b [nm] = ";
#      "RIN_f [s] = ";
#      "RIN_b [s] = "]


#      F.Valores = [I_980*1e3;
#           P0_980*1e3;
#      λ_0p_980*1e9;
#      Dλ_980*1e9;
#      I_1480*1e3;
#      P0_1480*1e3;
#      λ_0p_1480*1e9;
#      Dλ_1480*1e9;     
#      RL;
#      R0;
#      Z0;
#      Rho0;
#      T;
#      L;
#      1e3*Potf
#      1e3*Potb;
#      λmf*1e9;
#      λmb*1e9;
#      DλEff_f*1e9;
#      DλEff_b*1e9;
#      RIN_f;
#      RIN_f]

#      G=DataFrame()
#      G.L=vec(λ)
#      G.If=vec(Pf[:,Nz+1])
#      G.Ib=vec(Pb[:,1])



# #prPln("Salvar arquivo?  S/N")
#      tecla="N"
# if tecla == "S" begin
#      File=Dir*"/Resultados/M5-980B-L=8_5-T-300/"
#      File=File*string(now())

#      CSV.write(File*".csv", G)
#      CSV.write(File*".txt",F) 
#      savefig(File*".png")
#      end
# end

# println(F)
# println(I_980)
# ######End: Escrever os resultados calculados na tela e em arquivo ######

#**************************

"""
	end include Gorjian
"""
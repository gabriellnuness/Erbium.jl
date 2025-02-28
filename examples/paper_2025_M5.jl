using Erbium
using DelimitedFiles
using Interpolations



total_power_forward = zeros(7)
total_power_backward = zeros(7)
FC = 1

# FileExp_0="D:/Users/User/MyDrive/ITA/Data/3-Erbium_doped_fibers/2024_04_29_erbium_backward_forward/data/"
FileExp_0="S:/Users/Stinky/GD Academic/ITA/Data/3-Erbium_doped_fibers/2024_04_29_erbium_backward_forward/data/"


""" 980nm:  100 mA """ 

FileExp_B=FileExp_0*"experiment_3_m5_980nm_backward_100mA_1_99dBm.csv"
FileExp_F=FileExp_0*"experiment_4_m5_980nm_forward_100mA_1_03dBm.csv"

I_980 = 100e-3 # Corrente do laser 1480 nm

include("fiber_data.jl")
include("model_setup_integration.jl")

total_power_forward[1] = (1-RL)*Potf*1e3/FC
total_power_backward[1] = (1-R0)*Potb*1e3/FC










""" 980nm:  150 mA """ 

I_980 = 150e-3 #Corrente do laser 1480 nm
#subplot(242)
# figure(2)
FileExp_B=FileExp_0*"experiment_3_m5_980nm_backward_150mA_4_03dBm.csv"
FileExp_F=FileExp_0*"experiment_4_m5_980nm_forward_150mA_3_16dBm.csv"



include("fiber_data.jl")
include("model_setup_integration.jl") 

total_power_forward[2]=(1-RL)*Potf*1e3/FC
total_power_backward[2]=(1-R0)*Potb*1e3/FC








""" 980nm:  200 mA """ 

I_980 = 200e-3 #Corrente do laser 980 nm
#subplot(243)
figure(3)
FileExp_B=FileExp_0*"experiment_3_m5_980nm_backward_200mA_5_1dBm.csv"
FileExp_F=FileExp_0*"experiment_4_m5_980nm_forward_200mA_4_25dBm.csv"
include("fiber_data.jl")
include("model_setup_integration.jl")
total_power_forward[3]=(1-RL)*Potf*1e3/FC
total_power_backward[3]=(1-R0)*Potb*1e3/FC



""" 980nm:  250 mA """ 


I_980 = 250e-3 #Corrente do laser 980 nminclude("Dados_M5.jl")
#subplot(244)
figure(4)
FileExp_B=FileExp_0*"experiment_3_m5_980nm_backward_250mA_5_78dBm.csv"
FileExp_F=FileExp_0*"experiment_4_m5_980nm_forward_250mA_4_93dBm.csv"
include("fiber_data.jl")
include("model_setup_integration.jl")
total_power_forward[4]=(1-RL)*Potf*1e3/FC
total_power_backward[4]=(1-R0)*Potb*1e3/FC


""" 980nm:  300 mA """ 

I_980 = 300e-3 #Corrente do laser 980 nm
#subplot(245)
figure(5)
FileExp_B=FileExp_0*"experiment_3_m5_980nm_backward_300mA_6_27dBm.csv"
FileExp_F=FileExp_0*"experiment_4_m5_980nm_forward_300mA_5_41dBm.csv"
include("fiber_data.jl")
include("model_setup_integration.jl")
total_power_forward[5]=(1-RL)*Potf*1e3/FC
total_power_backward[5]=(1-R0)*Potb*1e3/FC


""" 980nm:  400 mA """ 

I_980 = 400e-3 #Corrente do laser 980 nm
#subplot(246)
figure(6)
FileExp_B=FileExp_0*"experiment_3_m5_980nm_backward_400mA_6_92dBm.csv"
FileExp_F=FileExp_0*"experiment_4_m5_980nm_forward_400mA_6_05dBm.csv"
include("fiber_data.jl")
include("model_setup_integration.jl")
total_power_forward[6]=(1-RL)*Potf*1e3/FC
total_power_backward[6]=(1-R0)*Potb*1e3/FC


""" 980nm:  500 mA """ 

I_980 = 500e-3 #Corrente do laser 980 nm
#subplot(247)
figure(7)
FileExp_B=FileExp_0*"experiment_3_m5_980nm_backward_500mA_7_35dBm.csv"
FileExp_F=FileExp_0*"experiment_4_m5_980nm_forward_500mA_6_46dBm.csv"
include("fiber_data.jl")
include("model_setup_integration.jl")
total_power_forward[7]=(1-RL)*Potf*1e3/FC
total_power_backward[7]=(1-R0)*Potb*1e3/FC

IP=[100;
150;
200;
250;
300;
400;
500]

PEF =[0.741;
0.988;
1.141;
1.245;
1.324;
1.431;
1.506]

PEF=[1.26765186585785;
2.07014134879104;
2.66072505979881;
3.11171633710602;
3.47536161443206;
4.02717034325459;
4.42588372362627]

PEB=[0.9;
1.2;
1.38;
1.514;
1.614;
1.754;
1.854]

PEB=[1.58124803927038;
2.52929799644615;
3.23593656929628;
3.78442584717093;
4.23642966049541;
4.92039535681451;
5.43250331492433]



global VarB=0
global VarF=0
for m=1:7
    global VarB
    global VarF
    VarB=VarB+(PEB[m]-total_power_backward[m])^2
    VarF=VarF+(PEF[m]-total_power_forward[m])^2
end

figure()
	plot(IP,total_power_backward)
	plot(IP,total_power_forward)
	plot(IP, PEB, "o", color=custom_plot_colors[1])
	plot(IP, PEF, "+", color=custom_plot_colors[2])
	yscale("log")


println([VarB VarF])
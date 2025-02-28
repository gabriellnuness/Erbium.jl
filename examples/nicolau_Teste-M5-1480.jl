using CSV
using DataFrames
using Dates
using DelimitedFiles
using Interpolations
using NativeFileDialog
using PyPlot
using PyCall
pygui(:gtk3)
using PyPlot
pygui(true)


PTF=zeros(7)
PTB=zeros(7)
FC=1

FileExp_0="/home/nicolau/Documentos/"
FileExp_0=FileExp_0*"IEAv/Fonte_erbio/Resultados Gabriel/Erbium/2024_04_29_erbium_backward_forward/data/"

#figure(1,figsize=(14,3))

I_1480 = 150e-3 #Corrente do laser 1480 nm
#subplot(241)
figure(1)
FileExp_B=FileExp_0*"experiment_10_m5_1480nm_backward_150mA_neg15_24dBm.txt"
FileExp_F=FileExp_0*"experiment_10_m5_1480nm_forward_150mA_neg6_72dBm.txt"
include("Dados-M5-1480.jl")
include("EDFA-Gorjian-M5-1480.jl")
PTF[1]=(1-RL)*Potf*1e3/FC
PTB[1]=(1-R0)*Potb*1e3/FC

I_1480 = 200e-3 #Corrente do laser 1480 nm
#subplot(242)
figure(2)
FileExp_B=FileExp_0*"experiment_10_m5_1480nm_backward_200mA_neg13_4dBm.txt"
FileExp_F=FileExp_0*"experiment_10_m5_1480nm_forward_200mA_neg3_84dBm.txt"
include("Dados-M5-1480.jl")
include("EDFA-Gorjian-M5-1480.jl")
PTF[2]=(1-RL)*Potf*1e3/FC
PTB[2]=(1-R0)*Potb*1e3/FC

I_1480 = 250e-3 #Corrente do laser 980 nm
#subplot(243)
figure(3)
FileExp_B=FileExp_0*"experiment_10_m5_1480nm_backward_250mA_neg12_35dBm.txt"
FileExp_F=FileExp_0*"experiment_10_m5_1480nm_forward_250mA_neg1_29dBm.txt"
1include("Dados-M5-1480.jl")
include("EDFA-Gorjian-M5-1480.jl")
PTF[3]=(1-RL)*Potf*1e3/FC
PTB[3]=(1-R0)*Potb*1e3/FC

I_1480 = 300e-3 #Corrente do laser 980 nminclude("Dados_M5.jl")
#subplot(244)
figure(4)
FileExp_B=FileExp_0*"experiment_10_m5_1480nm_backward_300mA_neg11_89dBm.txt"
FileExp_F=FileExp_0*"experiment_10_m5_1480nm_forward_300mA_0_2dBm.txt"
include("Dados-M5-1480.jl")
include("EDFA-Gorjian-M5-1480.jl")
PTF[4]=(1-RL)*Potf*1e3/FC
PTB[4]=(1-R0)*Potb*1e3/FC

I_1480 = 400e-3 #Corrente do laser 980 nm
#subplot(245)
figure(5)
FileExp_B=FileExp_0*"experiment_10_m5_1480nm_backward_400mA_neg11_46dBm.txt"
FileExp_F=FileExp_0*"experiment_10_m5_1480nm_forward_400mA_2_23dBm.txt"
include("Dados-M5-1480.jl")
include("EDFA-Gorjian-M5-1480.jl")
PTF[5]=(1-RL)*Potf*1e3/FC
PTB[5]=(1-R0)*Potb*1e3/FC

I_1480 = 500e-3 #Corrente do laser 980 nm
#subplot(246)
figure(6)
FileExp_B=FileExp_0*"experiment_10_m5_1480nm_backward_500mA_neg11_23dBm.txt"
FileExp_F=FileExp_0*"experiment_10_m5_1480nm_forward_500mA_3_54dBm.txt"
include("Dados-M5-1480.jl")
include("EDFA-Gorjian-M5-1480.jl")
PTF[6]=(1-RL)*Potf*1e3/FC
PTB[6]=(1-R0)*Potb*1e3/FC

I_1480 = 1000e-3 #Corrente do laser 980 nm
#subplot(247)
figure(7)
FileExp_B=FileExp_0*"experiment_10_m5_1480nm_backward_1000mA_neg11_02dBm.txt"
FileExp_F=FileExp_0*"experiment_10_m5_1480nm_forward_1000mA_6_99dBm.txt"
include("Dados-M5-1480.jl")
include("EDFA-Gorjian-M5-1480.jl")
PTF[7]=(1-RL)*Potf*1e3/FC
PTB[7]=(1-R0)*Potb*1e3/FC

IP=[150;
200;
250;
300;
400;
500;
1000]

PEF=[0.212813904598271;
0.413047501990161;
0.743019137896701;
1.0471285480509;
1.67109061431071;
2.25943577022098;
5.00034534976979]

PEB=[0.0299226463660819;
0.0457088189614875;
0.0582103217770872;
0.0647142615748583;
0.0714496326075513;
0.0753355563733717;
0.0790678627999825]



global VarB=0
global VarF=0
for m=1:7
    global VarB
    global VarF
    VarB=VarB+(PEB[m]-PTB[m])^2
    VarF=VarF+(PEF[m]-PTF[m])^2
end

#subplot(248)
figure(8)
plot(IP,PTB,"-b")
plot(IP,PTF,"-r")
plot(IP, PEB, "ob")
plot(IP, PEF, "+r")
yscale("log")
#grid("on")
println()
println()
println([VarB VarF])
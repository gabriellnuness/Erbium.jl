# using Pkg
# Pkg.add("https://github.com/gabriellnuness/Erbium.jl")
using Erbium

using DelimitedFiles
using Interpolations
using PyPlot
using DataFrames
using CSV


# sizes for Applied Optics paper
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["figure.autolayout"] = false  # Critical to get true size pdf
rcParams["savefig.bbox"] = "tight"
rcParams["savefig.pad_inches"] = 0
rcParams["font.size"] = 8
rcParams["axes.labelsize"] = 8
rcParams["legend.fontsize"] = 7
rcParams["axes.titlesize"] = 9
rcParams["figure.figsize"] = (8/2.54, 3.5/2.54)
rcParams["figure.dpi"] = 120
rcParams["xtick.major.size"] = 2.5
rcParams["xtick.minor.size"] = 1
rcParams["ytick.major.size"] = 2.5
rcParams["ytick.minor.size"] = 1
rcParams["xtick.major.width"] = 0.5
rcParams["xtick.minor.width"] = 0.3
rcParams["ytick.major.width"] = 0.5
rcParams["ytick.minor.width"] = 0.3
rcParams["axes.linewidth"] = 0.5


custom_plot_colors =   ["#6bc7f2"
                        "#fd5679"
                        "#f3bc6b"  
                        "#6cd776"
                        "#dee0de"]

PyPlot.pygui(true)




marker_step = 8
p1 = 101
p2 = 40

close("all")
"""
pumping 980
"""

PTF=zeros(7)
PTB=zeros(7)
FC=1
include("parameters.jl")

I_980 =  100e-3#200e-3 #Corrente do laser 1480 nm
I_1480 = 0.0 
L = 8.5 # metros - comprimento da fibra dopada com érbio ###

include("Dados-M5-980-A.jl")
include("EDFA-Gorjian-M5-980-A.jl")
PotCalcB = PotCalcB[p1:end-p2]
PotCalcF = PotCalcF[p1:end-p2]
PTF[1]=Potf*1e3/FC  #(1-RL)*Potf*1e3/FC
println(PTF)
PTB[1]=Potb*1e3/FC  #(1-R0)*Potb*1e3/FC
println(PTB)

figure(1)
plot(λ[p1:end-p2]*1e9, PotCalcB/MaxPotCalcB, color=custom_plot_colors[1], zorder=3)
plot(λ[p1:end-p2]*1e9, PotCalcF/MaxPotCalcF, color=custom_plot_colors[1], zorder=3) 
scatter(λ[p1:marker_step:end-p2].*1e9, PotCalcF[1:marker_step:end]./MaxPotCalcF, color=custom_plot_colors[1], marker="o", s=4) 
yscale("log")
ylim(1e-5,1.1)
xlim(1500, 1580)

# figure(2)
# plot(λ[p1:end-p2]*1e9, PotCalcF/MaxPotCalcF, color=custom_plot_colors[2], label="980 nm forward") 
# plot(λ[p1:end-p2]*1e9, PotCalcB/MaxPotCalcB, color=custom_plot_colors[1], label="backw980 nm backward")
# xlim(1500, 1580)
# ylim(0, 1.1)

# figure(3)
# plot(z,n2, color=custom_plot_colors[1], label="980")
# ylim(0, 1.1)

"""
pumping 1480
"""

PTF=zeros(7)
PTB=zeros(7)
FC=1
include("parameters.jl")

I_980 = 0.0 #500e-3 #Corrente do laser 1480 nm
I_1480 = 1.0 #3.0
L = 8.5 # metros - comprimento da fibra dopada com érbio ###

include("Dados-M5-980-A.jl")
include("EDFA-Gorjian-M5-980-A.jl")
PotCalcB = PotCalcB[p1:end-p2]
PotCalcF = PotCalcF[p1:end-p2]
PTF[1]=Potf*1e3/FC  #(1-RL)*Potf*1e3/FC
println(PTF)
PTB[1]=Potb*1e3/FC  #(1-R0)*Potb*1e3/FC

figure(1)
plot(λ[p1:end-p2]*1e9, PotCalcB/MaxPotCalcB, color=custom_plot_colors[3], "-", zorder=1)
plot(λ[p1:end-p2]*1e9, PotCalcF/MaxPotCalcF, color=custom_plot_colors[3], "-", zorder=1) 
scatter(λ[p1:marker_step:end-p2]*1e9, PotCalcF[1:marker_step:end]/MaxPotCalcF, color=custom_plot_colors[3], marker="D",s=2, zorder=1) 
yscale("log")
xlim(1500, 1580)

# figure(2)
# plot(λ[p1:end-p2]*1e9, PotCalcF/MaxPotCalcF, color=custom_plot_colors[2], "--", label="forward 1480") 
# plot(λ[p1:end-p2]*1e9, PotCalcB/MaxPotCalcB, color=custom_plot_colors[1], "--", label="backward 1480")
# xlim(1500, 1580)

# figure(3)
# plot(z,n2, color=custom_plot_colors[1], "--", label="1480")
# ylim(0, 1.1)

"""
pumping 980 +1480
"""
PTF=zeros(7)
PTB=zeros(7)
FC=1

I_980 =  100e-3#200e-3#500e-3 #Corrente do laser 1480 nm
I_1480 = 1.0
L = 8.5 # metros - comprimento da fibra dopada com érbio ###

include("Dados-M5-980-A.jl")
include("EDFA-Gorjian-M5-980-A.jl")
PotCalcB = PotCalcB[p1:end-p2]
PotCalcF = PotCalcF[p1:end-p2]
PTF[1]=Potf*1e3/FC  #(1-RL)*Potf*1e3/FC
PTB[1]=Potb*1e3/FC  #(1-R0)*Potb*1e3/FC

figure(1)
plot(λ[p1:end-p2]*1e9, PotCalcB/MaxPotCalcB, color=custom_plot_colors[2], "-", zorder=2)
plot(λ[p1:end-p2]*1e9, PotCalcF/MaxPotCalcF, color=custom_plot_colors[2], "-", zorder=2) 
scatter(λ[p1:marker_step:end-p2]*1e9, PotCalcF[1:marker_step:end]/MaxPotCalcF, color=custom_plot_colors[2], marker="s", s=4, zorder=2) 
yscale("log")
# xlim(1500, 1580)





# figure(2)
# plot(λ[p1:end-p2]*1e9, PotCalcF/MaxPotCalcF, color=custom_plot_colors[2], ":", label="forward dual") 
# plot(λ[p1:end-p2]*1e9, PotCalcB/MaxPotCalcB, color=custom_plot_colors[1], ":", label="backward dual")
# xlim(1500, 1580)

# figure(3)
# plot(z, n2, color=custom_plot_colors[1], ":", label="dual")
# ylim(0, 1.1)


figure(1)
xlabel("wavelength [nm]")
ylabel("normalized power")
plot([], [], color=custom_plot_colors[1], linestyle="-", label="980 nm backward")
plot([], [], color=custom_plot_colors[1], linestyle="-", marker="o", markersize=2, label="980 nm forward")
plot([], [], color=custom_plot_colors[2], linestyle="-", label="dual pump backward")
plot([], [], color=custom_plot_colors[2], linestyle="-", marker="s", markersize=2, label="dual pump forward")
plot([], [], color=custom_plot_colors[3], linestyle="-", label="1480 nm backward")
plot([], [], color=custom_plot_colors[3], linestyle="-", marker="D", markersize=2, label="1480 nm forward")
yscale("log")
ylim(1e-4, 1.1)
legend(frameon=false, loc="upper center", ncol=3, bbox_to_anchor=(0.5, -0.35), fontsize=6.5, handlelength=0.7, labelspacing=0.1, columnspacing=0.2,handletextpad=0.3)
savefig("paper_plots/output_figures/fig14_log.pdf", bbox_inches="tight", transparent=true, pad_inches=0.02)


# figure(2)
# legend()
# xlabel("wavelength [nm]")
# ylabel("normalized power")
# savefig("paper_plots/output_figures/fig14_lin.pdf", bbox_inches="tight", transparent=true, pad_inches=0)

# figure(3)
# legend()
# xlabel("fiber length [m]")
# ylabel("level 2 population")
# savefig("paper_plots/output_figures/fig14_population.pdf", bbox_inches="tight", transparent=true, pad_inches=0)












"""
generating figure 13 from paper
"""

"""
pumping 980
980 250 mA
"""
FC=1
L = 8.5
currents_980 = [0;50;100;150;200;250;300;400;500]*1e-3
PTF_980=zeros(size(currents_980))
PTB_980=zeros(size(currents_980))
mean_wavelength_980_b=zeros(size(currents_980))
mean_wavelength_980_f=zeros(size(currents_980))
bw_980_f=zeros(size(currents_980))
bw_980_b=zeros(size(currents_980))
I_980 = 0
I_1480 = 0

for (index, current) in enumerate(currents_980)

	I_980 =  current
	I_1480 = 1000e-3
	include("Dados-M5-980-A.jl")
	include("EDFA-Gorjian-M5-980-A.jl")
	PTF_980[index]=Potf*1e3/FC  #(1-RL)*Potf*1e3/FC
	PTB_980[index]=Potb*1e3/FC  #(1-R0)*Potb*1e3/FC
	mean_wavelength_980_f[index]=λmf
	mean_wavelength_980_b[index]=λmb
	bw_980_f[index]=DλEff_f
	bw_980_b[index]=DλEff_b
end

figure()
plot(currents_980, PTB_980)
plot(currents_980, PTF_980)
yscale("log")


figure()
plot(currents_980, mean_wavelength_980_b*1e9)
plot(currents_980, mean_wavelength_980_f*1e9)


figure()
plot(currents_980, bw_980_b*1e9)
plot(currents_980, bw_980_f*1e9)


"""
pumping 1480
"""

L = 8.5
currents_1480 = [0;150;200;250;300;400;500;700;1000]*1e-3
PTF_1480=zeros(size(currents_1480))
PTB_1480=zeros(size(currents_1480))
mean_wavelength_1480_b=zeros(size(currents_1480))
mean_wavelength_1480_f=zeros(size(currents_1480))
bw_1480_f=zeros(size(currents_1480))
bw_1480_b=zeros(size(currents_1480))
I_980 = 0
I_1480 = 0


for (index, current) in enumerate(currents_1480)

	I_980 =  250e-3
	I_1480 = current
	include("Dados-M5-980-A.jl")
	include("EDFA-Gorjian-M5-980-A.jl")
	PTF_1480[index]=Potf*1e3/FC  #(1-RL)*Potf*1e3/FC
	PTB_1480[index]=Potb*1e3/FC  #(1-R0)*Potb*1e3/FC
	mean_wavelength_1480_f[index]=λmf
	mean_wavelength_1480_b[index]=λmb
	bw_1480_f[index]=DλEff_f
	bw_1480_b[index]=DλEff_b
end

figure()
plot(currents_1480, PTB_1480)
plot(currents_1480, PTF_1480)
yscale("log")


figure()
plot(currents_1480, mean_wavelength_1480_b*1e9)
plot(currents_1480, mean_wavelength_1480_f*1e9)


figure()
plot(currents_1480, bw_1480_b*1e9)
plot(currents_1480, bw_1480_f*1e9)



df = DataFrame(	Current_980=currents_980,
				Power_980_forward=PTF_980,
				Power_980_backward=PTB_980,
				Bandwidth_980_forward=bw_980_f,
				Bandwidth_980_backward=bw_980_b,
				Mean_wavelength_980_forward=mean_wavelength_980_f,
				Mean_wavelength_980_backward=mean_wavelength_980_b,
				Current_1480=currents_1480,
				Power_1480_forward=PTF_1480,
				Power_1480_backward=PTB_1480,
				Bandwidth_1480_forward=bw_1480_f,
				Bandwidth_1480_backward=bw_1480_b,
				Mean_wavelength_1480_forward=mean_wavelength_1480_f,
				Mean_wavelength_1480_backward=mean_wavelength_1480_b,)
CSV.write("paper_plots/simulaion_dual_pumping_980_250ma.csv", df)







"""
generating figure 13 from paper  
980 100 mA
"""

"""
pumping 980
"""
FC=1
L = 8.5
currents_980 = [0;50;100;150;200;250;300;400;500]*1e-3
PTF_980=zeros(size(currents_980))
PTB_980=zeros(size(currents_980))
mean_wavelength_980_b=zeros(size(currents_980))
mean_wavelength_980_f=zeros(size(currents_980))
bw_980_f=zeros(size(currents_980))
bw_980_b=zeros(size(currents_980))
I_980 = 0
I_1480 = 0

for (index, current) in enumerate(currents_980)

	I_980 =  current
	I_1480 = 1000e-3
	include("Dados-M5-980-A.jl")
	include("EDFA-Gorjian-M5-980-A.jl")
	PTF_980[index]=Potf*1e3/FC  #(1-RL)*Potf*1e3/FC
	PTB_980[index]=Potb*1e3/FC  #(1-R0)*Potb*1e3/FC
	mean_wavelength_980_f[index]=λmf
	mean_wavelength_980_b[index]=λmb
	bw_980_f[index]=DλEff_f
	bw_980_b[index]=DλEff_b
end

figure()
plot(currents_980, PTB_980)
plot(currents_980, PTF_980)
yscale("log")


figure()
plot(currents_980, mean_wavelength_980_b*1e9)
plot(currents_980, mean_wavelength_980_f*1e9)


figure()
plot(currents_980, bw_980_b*1e9)
plot(currents_980, bw_980_f*1e9)




"""
pumping 1480
"""
L = 8.5
currents_1480 = [0;150;200;250;300;400;500;700;1000]*1e-3
PTF_1480=zeros(size(currents_1480))
PTB_1480=zeros(size(currents_1480))
mean_wavelength_1480_b=zeros(size(currents_1480))
mean_wavelength_1480_f=zeros(size(currents_1480))
bw_1480_f=zeros(size(currents_1480))
bw_1480_b=zeros(size(currents_1480))
I_980 = 0
I_1480 = 0


for (index, current) in enumerate(currents_1480)

	I_980 =  100e-3
	I_1480 = current
	include("Dados-M5-980-A.jl")
	include("EDFA-Gorjian-M5-980-A.jl")
	PTF_1480[index]=Potf*1e3/FC  #(1-RL)*Potf*1e3/FC
	PTB_1480[index]=Potb*1e3/FC  #(1-R0)*Potb*1e3/FC
	mean_wavelength_1480_f[index]=λmf
	mean_wavelength_1480_b[index]=λmb
	bw_1480_f[index]=DλEff_f
	bw_1480_b[index]=DλEff_b
end

figure()
plot(currents_1480, PTB_1480)
plot(currents_1480, PTF_1480)
yscale("log")


figure()
plot(currents_1480, mean_wavelength_1480_b*1e9)
plot(currents_1480, mean_wavelength_1480_f*1e9)


figure()
plot(currents_1480, bw_1480_b*1e9)
plot(currents_1480, bw_1480_f*1e9)



df = DataFrame(	Current_980=currents_980,
				Power_980_forward=PTF_980,
				Power_980_backward=PTB_980,
				Bandwidth_980_forward=bw_980_f,
				Bandwidth_980_backward=bw_980_b,
				Mean_wavelength_980_forward=mean_wavelength_980_f,
				Mean_wavelength_980_backward=mean_wavelength_980_b,
				Current_1480=currents_1480,
				Power_1480_forward=PTF_1480,
				Power_1480_backward=PTB_1480,
				Bandwidth_1480_forward=bw_1480_f,
				Bandwidth_1480_backward=bw_1480_b,
				Mean_wavelength_1480_forward=mean_wavelength_1480_f,
				Mean_wavelength_1480_backward=mean_wavelength_1480_b)
CSV.write("paper_plots/simulaion_dual_pumping_980_100ma.csv", df)







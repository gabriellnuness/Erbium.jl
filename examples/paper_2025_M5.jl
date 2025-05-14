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
	include("../src/fiber_data.jl")
	include("../src/core_integrator.jl")
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
	include("../src/fiber_data.jl")
	include("../src/core_integrator.jl")
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



df_100ma = DataFrame(	Current_980=currents_980,
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
CSV.write("simulaion_dual_pumping_980_100ma.csv", df)

benchmark_df = CSV.read("benchmark/simulaion_dual_pumping_980_100ma.csv", DataFrame)




"""
testing results
"""
# compare with benchmark values
function compare_with_tolerance(a, b; tol=1e-6)
    if eltype(a) <: AbstractFloat && eltype(b) <: AbstractFloat
        return all(abs.(a .- b) .< tol)
    else
        return a == b
    end
end
df = df_100ma
for col in names(df)
    if hasproperty(benchmark_df, col)
        if compare_with_tolerance(df[!, col], benchmark_df[!, col])
            println("Column $col matches")
        else
            println("Column $col differs")
            # You could add more detailed comparison here
        end
    else
        println("Column $col not found in benchmark")
    end
end






"""
generating figure 13 from paper  
980 250 mA
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
	include("../src/fiber_data.jl")
	include("../src/core_integrator.jl")
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
	include("../src/fiber_data.jl")
	include("../src/core_integrator.jl")
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



df_250ma = DataFrame(	Current_980=currents_980,
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
CSV.write("simulaion_dual_pumping_980_250ma.csv", df)

benchmark_df = CSV.read("benchmark/simulaion_dual_pumping_980_250ma.csv", DataFrame)

df = df_250ma
for col in names(df)
    if hasproperty(benchmark_df, col)
        if compare_with_tolerance(df[!, col], benchmark_df[!, col])
            println("Column $col matches")
        else
            println("Column $col differs")
            # You could add more detailed comparison here
        end
    else
        println("Column $col not found in benchmark")
    end
end


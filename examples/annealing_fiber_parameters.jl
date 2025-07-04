"""
Using Monte-Carlo simulation with annealing optimization 
to find the fiber parameters that are not exactly the calculated values 
from literature. For example, reflectivity on fiber ends, dopant level, temperature
"""

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



# mutable struct to input into integrator function and update in optimization algorithm
mutable struct FiberParameters
	RL_980::Float64
	RL_1480::Float64
	R0_980::Float64
	R0_1480::Float64
	T::Float64
	Ksigma::Float64
	Z::Float64
	ϵ1::Float64
	ϵ2::Float64
	gama::Float64
	G::Float64
	L::Float64
end

function FiberParameters(; # function with kwarg
# define some global initial values
	RL_980 = 0.0,
	RL_1480 = 0.3,
	R0_980 = 0.0,
	R0_1480 = 0.0,
	T = 273+27,
	Ksigma = 0.85,
	Z = 1.29e-9,
	ϵ1 =0.1,
	ϵ2 = 0.0159,
	gama = 0.01,
	G = 1.0,
	L = 8.5)
	
	# instantiate
	FiberParameters(RL_980, RL_1480, R0_980, R0_1480, T, Ksigma, Z, ϵ1, ϵ2, gama, G, L)
end

fiber_parameters = FiberParameters()



""" pumping 980 nm: 100 mA """

FC = 1
# currents_980 = [0;50;100;150;200;250;300;400;500]*1e-3
currents_980 = [500]*1e-3

# initialize variables for for-loop
PTF_980=zeros(size(currents_980))
PTB_980=zeros(size(currents_980))
mean_wavelength_980_b=zeros(size(currents_980))
mean_wavelength_980_f=zeros(size(currents_980))
bw_980_f=zeros(size(currents_980))
bw_980_b=zeros(size(currents_980))


@time begin


fiber_parameters = FiberParameters() # reboot to default values
	
figure()

# testing reflectivity change
for i=0:0.1:0.3
	println(i)

	# fiber_parameters.R0_980 = i
	# fiber_parameters.RL_980 = i
	# fiber_parameters.R0_1480 = i
	fiber_parameters.RL_1480 = i


	for (index, current) in enumerate(currents_980)


		println("\n\n sweeping 980: $(current*1000)mA\n\n")
		Nc = 75

		# this is a immutable tuple to input the data, for now it is fine because it is short,
		# but in future, transform in kwarg of integrate_single_pass() function
		input_parameters = (P0_980_F = pump_power_980(current),
							P0_1480_F = 0,
							P0_1480_B = pump_power_1480(1000e-3),
							P0_980_B = 0)

		# integrate spectrum
		spectra = integrator_single_pass(input_parameters,fiber_parameters)
		# characterize the generated spectrum
		params = calc_spec_parameters(spectra.power_forward[Nc:end], spectra.power_backward[Nc:end], spectra.λ[Nc:end])

		plot(spectra.power_backward[Nc:end],label="backward")
		plot(spectra.power_forward[Nc:end], label="forward")
		yscale("log")

		PTF_980[index]=params.Potf*1e3/FC  #(1-RL)*Potf*1e3/FC
		PTB_980[index]=params.Potb*1e3/FC  #(1-R0)*Potb*1e3/FC
		mean_wavelength_980_f[index]=params.λmf
		mean_wavelength_980_b[index]=params.λmb
		bw_980_f[index]=params.DλEff_f
		bw_980_b[index]=params.DλEff_b

	end # current loop - not in use for now, it was already implemented and I left it there...
end # reflectivity loop
end # @time

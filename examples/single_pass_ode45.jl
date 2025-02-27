using Erbium
using PyPlot
using OrdinaryDiffEq


@time begin
# Constants  
h = 6.6e-34
c = 299792458
Γ = 0.66   # optical fiber fill factor
fiber_loss = 0.0014    #  optical fiber loss in 1/m  

# Input data
λ_intervals = 100
z_intervals = 200
iter_intervals = 200
tol = 1e-3
λ1 = 1450e-9
λ2 = 1620e-9
dλ = (λ2-λ1)/λ_intervals

# Optical fiber info
fiber = "m5"
(β_abs0,β_emis0,τ21,τ3,fiber_diameter,NA,total_population,η) = optical_fiber(fiber)
temperature = 273+25
fiber_length = 10.0
dz = fiber_length / z_intervals


# separate spectrum in 3 parts
β_abs_pump = @view β_abs0[1:580,:] # view
β_emis_pump = zeros(size(β_abs_pump))
β_abs = @view β_abs0[581:end,:]
β_emis = @view β_emis0[581:end,:]





# Creating the pump spectrum
pump_power = 500e-3
ω0 = fiber_diameter/3
area = π*ω0^2
pump_intensity_peak = 2*pump_power/area
λc_pump = 974.4e-9
Δλ_pump = 2e-9
λ1_pump = 950e-9
λ2_pump = 1020e-9
dλ_pump = (λ2_pump-λ1_pump)/λ_intervals
λ_pump = range(start=λ1_pump, stop=λ2_pump, step=dλ_pump)
I_pump = generate_gaussian_spectrum(λ_pump, λc_pump, Δλ_pump, pump_intensity_peak)


# Initializing variables
λ_arr = range(start=λ1, stop=λ2, step=dλ)
z_arr = range(start=0, stop=z_intervals*dz, step=dz)
n2 = zeros(z_intervals+1)
I_forward = zeros(λ_intervals+1,z_intervals+1)
I_backward = zeros(λ_intervals+1,z_intervals+1)
I_pump_forward = zeros(λ_intervals+1,z_intervals+1)
I_pump_backward = zeros(λ_intervals+1,z_intervals+1)


# interpolate spectra to simulation data
β_abs_interp = linear_interpolation(β_abs[:,1], β_abs[:,2])
β_emis_interp = linear_interpolation(β_emis[:,1], β_emis[:,2])
β_abs_pump_interp = linear_interpolation(β_abs_pump[:,1], β_abs_pump[:,2])
β_emis_pump_interp = linear_interpolation(β_emis_pump[:,1], β_emis_pump[:,2])

# Defining spectral densities
convert_db_per_m_to_linear = log(10)/10
β_abs = convert_db_per_m_to_linear *  β_abs_interp.(λ_arr*1e9)
β_emis = convert_db_per_m_to_linear * β_emis_interp.(λ_arr*1e9)
β_abs_pump = convert_db_per_m_to_linear * β_abs_pump_interp.(λ_pump*1e9)
g = β_emis ./ trapz(λ_arr, β_emis) # ∫g(λ)⋅dλ = 1

# Pyplot commands
# figure()
# plot(λ_arr*1e9, β_abs)
# plot(λ_pump*1e9, β_abs_pump, color="tab:blue",label="_nolegend_")
# plot(λ_arr*1e9, β_emis)
# plot(λ_pump*1e9, I_pump./maximum(I_pump), color="tab:red")
#     legend(["Absorption", "Emission", "Pump"])


#####################
# Setup: Single pass
#####################
reflect_end = 0.0
reflect_begin = 1.0


# figure()
# plot(λ_pump,I_pump_forward[:,1])
# plot(λ_arr,I_forward[:,1])
#     legend(["Initial Pump","Initial Beam"])



# solving pump_forward with 5th order 
# the differential equation should be defined as f(u,t,p)
# u is the variable to solve,
# t is the time, or the integration base
# parameters are the constant terms for the ODE
function dIp_dz(I, parameters, z)
# TODO: include the n2 as a second function
# TODO: compare using just the first n2 to get the final Itensity with using all n2
    β_abs_pump, β_emis_pump, n2, Γ, fiber_loss = parameters
    pump_absorption = Γ*β_abs_pump*(1-n2)
    pump_emission =  Γ*β_emis_pump*n2
    
    return I*(-pump_absorption + pump_emission - fiber_loss)
end
function dI_dz(I, parameters, z)
    β_abs, β_emis, n2, Γ, fiber_loss, Rho = parameters
    absorption = Γ*β_abs*(1-n2)
    emission = Γ*β_emis*n2
    spontaneous = Rho*n2
    
    return I*(-absorption + emission - fiber_loss) + spontaneous
end


conv1 = 1
conv2 = 1
conv = ones(iter_intervals)

end

@time begin

for k = 1:iter_intervals

    """Forward integration --->"""
 
    # Initializing Intensitites for forward pass 
    I_pump_forward[:,1] = @. I_pump[:] + reflect_begin*I_pump_backward[:,1]
    I_forward[:,1] = @. reflect_begin * I_backward[:,1]

    @inbounds for λ in eachindex(λ_pump)
    # Find pump intensity for all λ from z=0 ->: z_end
    # and save all z steps in between
        parameters = (β_abs_pump[λ],β_emis_pump[λ],n2[1],Γ,fiber_loss) # prefer using tuple instead of arrays due to performance
        prob = ODEProblem(dIp_dz, I_pump_forward[λ,1], (0.0, fiber_length), parameters)
        sol = solve(prob, Tsit5(), saveat=0.05)

        I_pump_forward[λ,:] = sol.u
    end


    # the pump intensity is used to update the N2 level population
    for z = 1:z_intervals+1
        wp = 0
        wa = 0
        we = 0
        for λ = 1:λ_intervals+1
        # calculate W = ∫σ⋅I/(hν)dλ  --->  W = (1/NT) ⋅ ∫β⋅I⋅λ/(hc)dλ  --->   W = (1/NT) ⋅ ∫β⋅(If+Ib)⋅λ/(hc)dλ
            wp = wp + Γ*(I_pump_forward[λ,z]+I_pump_backward[λ,z])*λ_pump[λ]/(h*c)*β_abs_pump[λ]*(1-n2[z]) * dλ_pump/total_population
            wa = wa + Γ*(I_forward[λ,z]+I_backward[λ,z])*λ_arr[λ]/(h*c)*β_abs[λ]*(1-n2[z]) * dλ/total_population
            we = we + Γ*(I_forward[λ,z]+I_backward[λ,z])*λ_arr[λ]/(h*c)*β_emis[λ]*(1-n2[z]) * dλ/total_population
        end
        n2[z] = (wp+wa)*τ21/((wp+wa+we)*τ21+1)
    end

    # having the N2 population updated, the beam intensity is calculated
    @inbounds for λ in eachindex(λ_arr)
    # calculate beam intensity for all λ from z=0 ->: z_end
        Rho = h*c/λ_arr[λ]*g[λ]*η*0.5*(1-(1-NA^2)^(1/2))*total_population/τ21
        parameters = (β_abs[λ],β_emis[λ],n2[1],Γ,fiber_loss, Rho) # prefer using tuple instead of arrays due to performance
        prob = ODEProblem(dI_dz, I_forward[λ,1], (0.0, fiber_length), parameters)
        sol = solve(prob, Tsit5(), saveat=0.05)

        I_forward[λ,:] = sol.u
    end
    
    
    """<---- Backward"""
    # initializing intensities for backward pass
    I_pump_backward[:,z_intervals+1] = @. reflect_end * I_pump_forward[:,z_intervals+1]
    I_backward[:,z_intervals+1] = @. reflect_end * I_forward[:,z_intervals+1]

    @inbounds for λ in eachindex(λ_pump)
    # Find pump intensity for all λ from z=0 ->: z_end
    # and save all z steps in between
        parameters = (β_abs_pump[λ],β_emis_pump[λ],n2[1],Γ,fiber_loss) # prefer using tuple instead of arrays due to performance
        prob = ODEProblem(dIp_dz, I_pump_backward[λ,1], (0.0, fiber_length), parameters)
        sol = solve(prob, Tsit5(), saveat=0.05)

        I_pump_backward[λ,:] = sol.u
    end

    # update all N2 population after backward pass
    for j=1:z_intervals+1
        wp = 0
        wa = 0
        we = 0
        for i=1:λ_intervals+1
            wp = wp + Γ*(I_pump_forward[i,j]+I_pump_backward[i,j])*λ_pump[i]/(h*c)*β_abs_pump[i]*(1 - n2[j])*dλ_pump/total_population
            wa = wa + Γ*(I_forward[i,j]+I_backward[i,j])*λ_arr[i]/(h*c)*β_abs[i]*(1 - n2[j])*dλ/total_population
            we = we + Γ*(I_forward[i,j]+I_backward[i,j])*λ_arr[i]/(h*c)*β_emis[i]*(1 - n2[j])*dλ/total_population
        end
        n2[j] = (wp+wa)*τ21 / ((wp+wa+we)*τ21+1)
    end


    # calculate the backwar beam intensity after the n2 population was updated
    # @inbounds for λ in eachindex(λ_arr)
    # # calculate beam intensity for all λ from z=0 ->: z_end
    #     Rho = h*c/λ_arr[λ]*g[λ]*η*0.5*(1-(1-NA^2)^(1/2))*total_population/τ21
    #     parameters = (β_abs[λ],β_emis[λ],n2[end],Γ,fiber_loss, Rho) # prefer using tuple instead of arrays due to performance
    #     prob = ODEProblem(dI_dz, I_backward[λ,1], (0.0,fiber_length), parameters)
    #     sol = solve(prob, Tsit5(), saveat=0.05)

    #     I_backward[λ,:] = sol.u
    # end
     for j = 1:z_intervals    
            jj=z_intervals+1-j
            for i=1:λ_intervals+1
                # rdIdz(Int,Γ,β_abs,n2,β_emis,fiber_loss,rho,τ21)
                Rho=h*c/λ_arr[i]*g[i]*η*0.5*(1-(1-NA^2)^(1/2))*total_population/τ21
                fp1 = dIdz(I_backward[i,jj], Γ,β_abs[i],n2[jj],β_emis[i],fiber_loss,Rho) * dz
                fp2 = dIdz(I_backward[i,jj]+fp1/2, Γ,β_abs[i],n2[jj],β_emis[i],fiber_loss,Rho) * dz
                fp3 = dIdz(I_backward[i,jj]+fp2/2, Γ,β_abs[i],n2[jj],β_emis[i],fiber_loss,Rho) * dz
                fp4 = dIdz(I_backward[i,jj]+fp3, Γ,β_abs[i],n2[jj],β_emis[i],fiber_loss,Rho) * dz
                I_backward[i,jj]=I_backward[i,jj+1]+(fp1 + 2*fp2 + 2*fp3 + fp4)/6
            end
    end   


    # Checking convergence
    println("iteration: $k")
    conv2 = maximum(I_forward)
    conv[k+1] = (abs((conv2-conv1)/conv1) )
    if conv[k+1] > 10
        conv[k+1] = 10
    end
    println("convergence: $(conv[k+1])")
    if k > 3
        if conv[k-1]<tol && conv[k]<tol && conv[k+1]<tol
            break
        end
    end     
    conv1 = conv2

end

end # time

figure()
plot(λ_arr*1e9, I_forward[:,end])
plot(λ_arr*1e9, I_backward[:,1])
    xlabel("Wavelength [nm]")
    ylabel("Intensity [W/m²]")
    yscale("log")

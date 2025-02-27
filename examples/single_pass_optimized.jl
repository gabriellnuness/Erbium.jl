using Erbium
using PyPlot

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
fiber_length = 10
dz = fiber_length/z_intervals
reflect_end = 0.0
reflect_begin = 1.0

# separate spectrum in 3 parts
β_abs_pump = β_abs0[1:580,:]
β_emis_pump = zeros(size(β_abs_pump))
β_abs = β_abs0[581:end,:]
β_emis = β_emis0[581:end,:]





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
λarr = range(start=λ1, stop=λ2, step=dλ)
zarr = range(start=0, stop=z_intervals*dz, step=dz)
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
β_abs = convert_db_per_m_to_linear *  β_abs_interp.(λarr*1e9)
β_emis = convert_db_per_m_to_linear * β_emis_interp.(λarr*1e9)
β_abs_pump = convert_db_per_m_to_linear * β_abs_pump_interp.(λ_pump*1e9)
g = β_emis ./ trapz(λarr, β_emis) # ∫g(λarr)⋅dλ = 1


figure()
plot(λarr*1e9, β_abs)
plot(λ_pump*1e9, β_abs_pump, color="tab:blue",label="_nolegend_")
plot(λarr*1e9, β_emis)
plot(λ_pump*1e9, I_pump./maximum(I_pump), color="tab:red")
    legend(["Absorption", "Emission", "Pump"])






conv1 = 1
conv2 = 1
conv = ones(iter_intervals)
end


@time begin
for k = 1:iter_intervals

    # direction --->
    # initial conditions
    for λ = 1:λ_intervals+1
        I_pump_forward[λ,1] = I_pump[λ] + reflect_begin*I_pump_backward[λ,1]
        I_forward[λ,1] = reflect_begin * I_backward[λ,1]
    end 

    for z = 1:z_intervals
        for λ = 1:λ_intervals+1
        # runge kutta!!!
            # res=dIpdz_opt(Int,Γ,β_abs_pump,n2,β_emis_pump,fiber_loss)
            fp1 = dIpdz_opt(I_pump_forward[λ,z],n2[z],Γ,β_abs_pump[λ],β_emis_pump[λ],fiber_loss)*dz
            fp2 = dIpdz_opt(I_pump_forward[λ,z]+fp1/2,n2[z],Γ,β_abs_pump[λ],β_emis_pump[λ],fiber_loss)*dz
            fp3 = dIpdz_opt(I_pump_forward[λ,z]+fp2/2,n2[z],Γ,β_abs_pump[λ],β_emis_pump[λ],fiber_loss)*dz          
            fp4 = dIpdz_opt(I_pump_forward[λ,z]+fp3,n2[z],Γ,β_abs_pump[λ],β_emis_pump[λ],fiber_loss)*dz                    
            I_pump_forward[λ,z+1] = I_pump_forward[λ,z] + ( fp1 + 2*fp2 + 2*fp3 + fp4)/6
        end
    end


    for z = 1:z_intervals+1
    # calculate n2 for every zarr
        wp = 0
        wa = 0
        we = 0
        for λ = 1:λ_intervals+1
        # calculate W = ∫σ⋅I/(hν)dλ 
        # W = (1/NT) ⋅ ∫β⋅I⋅λ/(hc)dλ 
        # W = (1/NT) ⋅ ∫β⋅(If+Ib)⋅λ/(hc)dλ
            wp = wp + Γ*(I_pump_forward[λ,z]+I_pump_backward[λ,z])*λ_pump[λ]/(h*c)*β_abs_pump[λ]*(1-n2[z]) * dλ_pump/total_population
            wa = wa + Γ*(I_forward[λ,z]+I_backward[λ,z])*λarr[λ]/(h*c)*β_abs[λ]*(1-n2[z]) * dλ/total_population
            we = we + Γ*(I_forward[λ,z]+I_backward[λ,z])*λarr[λ]/(h*c)*β_emis[λ]*(1-n2[z]) * dλ/total_population
        end
        n2[z] = (wp+wa)*τ21/((wp+wa+we)*τ21+1)
    end

    for zj = 1:z_intervals          
            for λi = 1:λ_intervals+1
            # calculate ∫dI(λarr,zarr)dz
                Rho = h*c/λarr[λi]*g[λi]*η*0.5*(1-(1-NA^2)^(1/2))*total_population/τ21
                p = (Γ,β_abs[λi],β_emis[λi],fiber_loss,Rho)
                fp1 = dIdz_opt(I_forward[λi,zj],        n2[zj], p)*dz
                fp2 = dIdz_opt(I_forward[λi,zj]+fp1/2,  n2[zj], p)*dz
                fp3 = dIdz_opt(I_forward[λi,zj]+fp2/2,  n2[zj], p)*dz
                fp4 = dIdz_opt(I_forward[λi,zj]+fp3,    n2[zj], p)*dz
                I_forward[λi,zj+1] = I_forward[λi,zj]+(fp1+ 2*fp2 + 2*fp3 + fp4)/6
            end
    end 

    
    # <---- direction
    # initial conditions
    for λ = 1:λ_intervals+1
        I_pump_backward[λ,z_intervals+1] = reflect_end * I_pump_forward[λ,z_intervals+1]
        I_backward[λ,z_intervals+1] = reflect_end * I_forward[λ,z_intervals+1]
    end #  λ

    for z = 1:z_intervals
        for λ = 1:λ_intervals+1
            # res=dIpdz_opt(Int,Γ,β_abs_pump,n2,fiber_loss)
            fp1= dIpdz_opt(I_pump_forward[λ,z],n2[z],Γ,β_abs_pump[λ],β_emis_pump[λ],fiber_loss)*dz
            fp2= dIpdz_opt(I_pump_forward[λ,z]+fp1/2,n2[z],Γ,β_abs_pump[λ],β_emis_pump[λ],fiber_loss)*dz
            fp3= dIpdz_opt(I_pump_forward[λ,z]+fp2/2,n2[z],Γ,β_abs_pump[λ],β_emis_pump[λ],fiber_loss)*dz          
            fp4= dIpdz_opt(I_pump_forward[λ,z]+fp3,n2[z],Γ,β_abs_pump[λ],β_emis_pump[λ],fiber_loss)*dz                    
            I_pump_forward[λ,z+1] = I_pump_forward[λ,z] + ( fp1 + 2*fp2 + 2*fp3 + fp4)/6
        end
    end

    for z=1:z_intervals+1
        wp = 0
        wa = 0
        we = 0
        for λ=1:λ_intervals+1
            wp = wp + Γ*(I_pump_forward[λ,z]+I_pump_backward[λ,z])*λ_pump[λ]/(h*c)*β_abs_pump[λ]*(1 - n2[z])*dλ_pump/total_population
            wa = wa + Γ*(I_forward[λ,z]+I_backward[λ,z])*λarr[λ]/(h*c)*β_abs[λ]*(1 - n2[z])*dλ/total_population
            we = we + Γ*(I_forward[λ,z]+I_backward[λ,z])*λarr[λ]/(h*c)*β_emis[λ]*(1 - n2[z])*dλ/total_population
        end
        n2[z] = (wp+wa)*τ21 / ((wp+wa+we)*τ21+1)
    end

    for z = 1:z_intervals    
            z_back=z_intervals+1-z
            for λ=1:λ_intervals+1
                # rdIdz_opt(Int,Γ,β_abs,n2,β_emis,fiber_loss,rho,τ21)
                Rho=h*c/λarr[λ]*g[λ]*η*0.5*(1-(1-NA^2)^(1/2))*total_population/τ21
                p = (Γ,β_abs[λ],β_emis[λ],fiber_loss,Rho)
                fp1 = dIdz_opt(I_backward[λ,z_back],        n2[z_back],p) * dz
                fp2 = dIdz_opt(I_backward[λ,z_back]+fp1/2,  n2[z_back],p) * dz
                fp3 = dIdz_opt(I_backward[λ,z_back]+fp2/2,  n2[z_back],p) * dz
                fp4 = dIdz_opt(I_backward[λ,z_back]+fp3,    n2[z_back],p) * dz
                I_backward[λ,z_back]=I_backward[λ,z_back+1]+(fp1 + 2*fp2 + 2*fp3 + fp4)/6
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



figure()
plot(λarr*1e9, I_forward[:,end],alpha=0.2,linewidth=3)
plot(λarr*1e9, I_backward[:,1],alpha=0.2,linewidth=3)
    xlabel("Wavelength [nm]")
    ylabel("Intensity [W/m²]")
    yscale("log")

end
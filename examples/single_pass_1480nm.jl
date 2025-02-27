using Erbium
using PyPlot




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

# separate spectrum
β_abs = @view β_abs0[581:end,:]
β_emis = @view β_emis0[581:end,:]





# Creating the pump spectrum
pump_power = 500e-3
ω0 = fiber_diameter/3
area = π*ω0^2
pump_intensity_peak = 2*pump_power/area
λc_pump = 1480e-9 # 974.4e-9
Δλ_pump = 2e-9
λ1_pump = 950e-9
λ2_pump = 1020e-9
dλ_pump = (λ2_pump-λ1_pump)/λ_intervals
λ_pump = range(start=λ1_pump, stop=λ2_pump, step=dλ_pump)
I_pump = generate_gaussian_spectrum(λ_pump, λc_pump, Δλ_pump, pump_intensity_peak)


# Initializing variables
λ = range(start=λ1, stop=λ2, step=dλ)
z = range(start=0, stop=z_intervals*dz, step=dz)
n2 = zeros(z_intervals+1)
I_forward = zeros(λ_intervals+1,z_intervals+1)
I_backward = zeros(λ_intervals+1,z_intervals+1)


# interpolate spectra to simulation data
β_abs_interp = linear_interpolation(β_abs[:,1], β_abs[:,2])
β_emis_interp = linear_interpolation(β_emis[:,1], β_emis[:,2])

# Defining spectral densities

convert_db_per_m_to_linear = log(10)/10
β_abs = convert_db_per_m_to_linear *  β_abs_interp.(λ*1e9)
β_emis = convert_db_per_m_to_linear * β_emis_interp.(λ*1e9)
β_abs_pump = convert_db_per_m_to_linear * β_abs_pump_interp.(λ_pump*1e9)
g = β_emis ./ trapz(λ, β_emis) # ∫g(λ)⋅dλ = 1


figure()
plot(λ*1e9, β_abs)
plot(λ_pump*1e9, β_abs_pump, color="tab:blue",label="_nolegend_")
plot(λ*1e9, β_emis)
plot(λ_pump*1e9, I_pump./maximum(I_pump), color="tab:red")
    legend(["Absorption", "Emission", "Pump"])



conv1 = 1
conv2 = 1
conv = ones(iter_intervals)

@time begin

for k = 1:iter_intervals

    
    # direction --->
    # initial conditions
    for i = 1:λ_intervals+1
        I_pump_forward[i,1] = I_pump[i] + reflect_begin*I_pump_backward[i,1]
        I_forward[i,1] = reflect_begin * I_backward[i,1]
    end 


    for j = 1:z_intervals
        for i = 1:λ_intervals+1
            fp1 = dIpdz(I_pump_forward[i,j],         Γ,β_abs_pump[i],n2[j],β_emis_pump[i],fiber_loss)*dz
            fp2 = dIpdz(I_pump_forward[i,j]+fp1/2,   Γ,β_abs_pump[i],n2[j],β_emis_pump[i],fiber_loss)*dz
            fp3 = dIpdz(I_pump_forward[i,j]+fp2/2,   Γ,β_abs_pump[i],n2[j],β_emis_pump[i],fiber_loss)*dz          
            fp4 = dIpdz(I_pump_forward[i,j]+fp3,     Γ,β_abs_pump[i],n2[j],β_emis_pump[i],fiber_loss)*dz                    
            I_pump_forward[i,j+1] = I_pump_forward[i,j] + ( fp1 + 2*fp2 + 2*fp3 + fp4)/6
        end
    end


    for j = 1:z_intervals+1
    # calculate n2 for every z
        wp = 0
        wa = 0
        we = 0
        for i = 1:λ_intervals+1
        # calculate W = ∫σ⋅I/(hν)dλ 
        # W = (1/NT) ⋅ ∫β⋅I⋅λ/(hc)dλ 
        # W = (1/NT) ⋅ ∫β⋅(If+Ib)⋅λ/(hc)dλ
            wp = wp + Γ*(I_pump_forward[i,j]+I_pump_backward[i,j])*λ_pump[i]/(h*c)*β_abs_pump[i]*(1-n2[j]) * dλ_pump/total_population
            wa = wa + Γ*(I_forward[i,j]+I_backward[i,j])*λ[i]/(h*c)*β_abs[i]*(1-n2[j]) * dλ/total_population
            we = we + Γ*(I_forward[i,j]+I_backward[i,j])*λ[i]/(h*c)*β_emis[i]*(1-n2[j]) * dλ/total_population
        end
        n2[j] = (wp+wa)*τ21/((wp+wa+we)*τ21+1)
    end

    for zj = 1:z_intervals          
            for i = 1:λ_intervals+1
            # calculate ∫dI(λ,z)dz
                Rho = h*c/λ[i]*g[i]*η*0.5*(1-(1-NA^2)^(1/2))*total_population/τ21
                fp1 = dIdz(I_forward[i,zj],          Γ,β_abs[i],n2[zj],β_emis[i],fiber_loss,Rho)*dz
                fp2 = dIdz(I_forward[i,zj]+fp1/2,    Γ,β_abs[i],n2[zj],β_emis[i],fiber_loss,Rho)*dz
                fp3 = dIdz(I_forward[i,zj]+fp2/2,    Γ,β_abs[i],n2[zj],β_emis[i],fiber_loss,Rho)*dz
                fp4 = dIdz(I_forward[i,zj]+fp3,      Γ,β_abs[i],n2[zj],β_emis[i],fiber_loss,Rho)*dz
                I_forward[i,zj+1] = I_forward[i,zj]+(fp1+ 2*fp2 + 2*fp3 + fp4)/6
            end
    end 

    
    # <---- direction
    # initial conditions
    for i = 1:λ_intervals+1
        I_pump_backward[i,z_intervals+1] = reflect_end * I_pump_forward[i,z_intervals+1]
        I_backward[i,z_intervals+1] = reflect_end * I_forward[i,z_intervals+1]
    end #  i

    for j = 1:z_intervals
        for i = 1:λ_intervals+1
            # res=dIpdz(Int,Γ,β_abs_pump,n2,fiber_loss)
            fp1= dIpdz(I_pump_forward[i,j],Γ,β_abs_pump[i],n2[j],β_emis_pump[i],fiber_loss)*dz
            fp2= dIpdz(I_pump_forward[i,j]+fp1/2,Γ,β_abs_pump[i],n2[j],β_emis_pump[i],fiber_loss)*dz
            fp3= dIpdz(I_pump_forward[i,j]+fp2/2,Γ,β_abs_pump[i],n2[j],β_emis_pump[i],fiber_loss)*dz          
            fp4= dIpdz(I_pump_forward[i,j]+fp3,Γ,β_abs_pump[i],n2[j],β_emis_pump[i],fiber_loss)*dz                    
            I_pump_forward[i,j+1] = I_pump_forward[i,j] + ( fp1 + 2*fp2 + 2*fp3 + fp4)/6
        end
    end

    for j=1:z_intervals+1
        wp = 0
        wa = 0
        we = 0
        for i=1:λ_intervals+1
            wp = wp + Γ*(I_pump_forward[i,j]+I_pump_backward[i,j])*λ_pump[i]/(h*c)*β_abs_pump[i]*(1 - n2[j])*dλ_pump/total_population
            wa = wa + Γ*(I_forward[i,j]+I_backward[i,j])*λ[i]/(h*c)*β_abs[i]*(1 - n2[j])*dλ/total_population
            we = we + Γ*(I_forward[i,j]+I_backward[i,j])*λ[i]/(h*c)*β_emis[i]*(1 - n2[j])*dλ/total_population
        end
        n2[j] = (wp+wa)*τ21 / ((wp+wa+we)*τ21+1)
    end

    for j = 1:z_intervals    
            jj=z_intervals+1-j
            for i=1:λ_intervals+1
                # rdIdz(Int,Γ,β_abs,n2,β_emis,fiber_loss,rho,τ21)
                Rho=h*c/λ[i]*g[i]*η*0.5*(1-(1-NA^2)^(1/2))*total_population/τ21
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

end #time


figure()
plot(λ*1e9, I_forward[:,end],alpha=0.2,linewidth=2)
plot(λ*1e9, I_backward[:,1],alpha=0.2,linewidth=2)
    xlabel("Wavelength [nm]")
    ylabel("Intensity [W/m²]")
    yscale("log")


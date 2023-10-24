using Erbium
using Test
using PyPlot
using DelimitedFiles
using Interpolations
using Trapz

cd("..")

@testset "Spectrum simulation 980nm pump" begin

    # Constantes  
    h = 6.6e-34;
    c = 299792458;

    # Input data
    λ_intervals = 100
    z_intervals = 200
    iter_intervals = 500
    tol = 1e-3
    λ1 = 1450e-9
    λ2 = 1620e-9
    dλ = (λ2-λ1)/λ_intervals
    




    # Optical fiber info
    fiber = "m5"
    (β_abs0,β_emis0,τ21,τ3,fiber_diameter,NA,total_population,η) = optical_fiber(fiber)
    temperature = 273+25
    length = 10
    dz = length/z_intervals
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
    λc_pump = 980e-9
    Δλ_pump = 5e-9
    λ1_pump = 950e-9
    λ2_pump = 1020e-9
    dλ_pump = (λ2_pump-λ1_pump)/λ_intervals


    # Initializing variables
    λ_pump = range(start=λ1_pump, stop=λ2_pump, step=dλ_pump)
    λ = range(start=λ1, stop=λ2, step=dλ)
    z = range(start=0, stop=z_intervals*dz, step=dz)
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
    Γ = 0.66   # optical fiber fill factor
    gamma = 0.0014    #  optical fiber loss in 1/m  

    
    convert_db_per_m_to_linear = log(10)/10
    β_abs = convert_db_per_m_to_linear *  β_abs_interp.(λ*1e9)
    β_emis = convert_db_per_m_to_linear * β_emis_interp.(λ*1e9)
    β_abs_pump = convert_db_per_m_to_linear * β_abs_pump_interp.(λ_pump*1e9)
    g = β_emis ./ trapz(λ, β_emis) # ∫g(λ)⋅dλ = 1

    I_pump = @.  √(2/π) * pump_intensity_peak / dλ_pump * 
        exp(-(((λ_pump-λc_pump) / Δλ_pump)^2)) # Gaussiana distribution

        figure()
        plot(β_abs0[:,1], β_abs0[:,2],"o")
        plot(β_emis0[:,1], β_emis0[:,2],"o")
        plot(λ*1e9, β_abs,".")
        plot(λ*1e9, β_emis,".")
        plot(λ_pump*1e9, β_abs_pump,".")
    
        






    """Differential equation solver"""
    dIdz(I,Γ,β_abs,n2,β_emis,gamma,Rho) = I*(-Γ*β_abs*(1-n2)+Γ*β_emis*n2-gamma)+Rho*n2;
    dIpdz(I,Γ,β_abs_pump,n2,β_emis_pump,gamma) = I*(-Γ*β_abs_pump*(1-n2)+Γ*β_emis_pump*n2-gamma);

    Refa = 1
    Refb = 1
    conv = ones(iter_intervals)
    for k = 1:iter_intervals

        #  direction --->
        # initial conditions
        for i = 1:λ_intervals+1
            I_pump_forward[i,1] = I_pump[i] + reflect_begin*I_pump_backward[i,1]
            I_forward[i,1] = reflect_begin * I_backward[i,1]
        end 

        for j = 1:z_intervals
            for i = 1:λ_intervals+1
            # runge kutta!!!
                # res=dIpdz(Int,Γ,β_abs_pump,n2,β_emis_pump,gamma)
                fp1= dIpdz(I_pump_forward[i,j],Γ,β_abs_pump[i],n2[j],β_emis_pump[i],gamma)*dz
                fp2= dIpdz(I_pump_forward[i,j]+fp1/2,Γ,β_abs_pump[i],n2[j],β_emis_pump[i],gamma)*dz
                fp3= dIpdz(I_pump_forward[i,j]+fp2/2,Γ,β_abs_pump[i],n2[j],β_emis_pump[i],gamma)*dz          
                fp4= dIpdz(I_pump_forward[i,j]+fp3,Γ,β_abs_pump[i],n2[j],β_emis_pump[i],gamma)*dz                    
                I_pump_forward[i,j+1] = I_pump_forward[i,j] + ( fp1 + 2*fp2 + 2*fp3 + fp4)/6
            end
        end
        
        for j = 1:z_intervals+1
            wp=0; wa=0; we=0
            for i=1:λ_intervals+1
                wp=wp+Γ*(I_pump_forward[i,j]+I_pump_backward[i,j])*λ_pump[i]/(h*c)*β_abs_pump[i]*(1-n2[j])*dλ_pump/total_population
                wa=wa+Γ*(I_forward[i,j]+I_backward[i,j])*λ[i]/(h*c)*β_abs[i]*(1-n2[j])*dλ/total_population
                we=we+Γ*(I_forward[i,j]+I_backward[i,j])*λ[i]/(h*c)*β_emis[i]*(1-n2[j])*dλ/total_population
            end # i
            n2[j]=(wp+wa)*τ21/((wp+wa+we)*τ21+1)
        end # j

        for j=1:z_intervals          
                for i=1:λ_intervals+1
                    Rho = h*c/λ[i]*g[i]*η*0.5*(1-(1-NA^2)^(1/2))*total_population/τ21
                    fp1 = dIdz(I_forward[i,j],Γ,β_abs[i],n2[j],β_emis[i],gamma,Rho)*dz
                    fp2 = dIdz(I_forward[i,j]+fp1/2,Γ,β_abs[i],n2[j],β_emis[i],gamma,Rho)*dz
                    fp3 = dIdz(I_forward[i,j]+fp2/2,Γ,β_abs[i],n2[j],β_emis[i],gamma,Rho)*dz
                    fp4 = dIdz(I_forward[i,j]+fp3,Γ,β_abs[i],n2[j],β_emis[i],gamma,Rho)*dz
                    I_forward[i,j+1]=I_forward[i,j]+(fp1+ 2*fp2 + 2*fp3 + fp4)/6
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
                # res=dIpdz(Int,Γ,β_abs_pump,n2,gamma)
                fp1= dIpdz(I_pump_forward[i,j],Γ,β_abs_pump[i],n2[j],β_emis_pump[i],gamma)*dz
                fp2= dIpdz(I_pump_forward[i,j]+fp1/2,Γ,β_abs_pump[i],n2[j],β_emis_pump[i],gamma)*dz
                fp3= dIpdz(I_pump_forward[i,j]+fp2/2,Γ,β_abs_pump[i],n2[j],β_emis_pump[i],gamma)*dz          
                fp4= dIpdz(I_pump_forward[i,j]+fp3,Γ,β_abs_pump[i],n2[j],β_emis_pump[i],gamma)*dz                    
                I_pump_forward[i,j+1] = I_pump_forward[i,j] + ( fp1 + 2*fp2 + 2*fp3 + fp4)/6
            end
        end
    
        for j=1:z_intervals+1
            wp=0; wa=0; we=0
            for i=1:λ_intervals+1
                wp = wp + Γ*(I_pump_forward[i,j]+I_pump_backward[i,j])*λ_pump[i]/(h*c)*β_abs_pump[i]*(1 - n2[j])*dλ_pump/total_population
                wa = wa + Γ*(I_forward[i,j]+I_backward[i,j])*λ[i]/(h*c)*β_abs[i]*(1 - n2[j])*dλ/total_population
                we = we + Γ*(I_forward[i,j]+I_backward[i,j])*λ[i]/(h*c)*β_emis[i]*(1 - n2[j])*dλ/total_population
            end # i
            n2[j] = (wp+wa)*τ21 / ((wp+wa+we)*τ21+1)
        end
    
        for j = 1:z_intervals    
                jj=z_intervals+1-j
                for i=1:λ_intervals+1
                    # rdIdz(Int,Γ,β_abs,n2,β_emis,gamma,rho,τ21)
                    Rho=h*c/λ[i]*g[i]*η*0.5*(1-(1-NA^2)^(1/2))*total_population/τ21
                    fp1 = dIdz(I_backward[i,jj], Γ,β_abs[i],n2[jj],β_emis[i],gamma,Rho) * dz
                    fp2 = dIdz(I_backward[i,jj]+fp1/2, Γ,β_abs[i],n2[jj],β_emis[i],gamma,Rho) * dz
                    fp3 = dIdz(I_backward[i,jj]+fp2/2, Γ,β_abs[i],n2[jj],β_emis[i],gamma,Rho) * dz
                    fp4 = dIdz(I_backward[i,jj]+fp3, Γ,β_abs[i],n2[jj],β_emis[i],gamma,Rho) * dz
                    I_backward[i,jj]=I_backward[i,jj+1]+(fp1 + 2*fp2 + 2*fp3 + fp4)/6
                end
        end # j     



        # Checking convergence
        println("iteration: $k")
        Refb = maximum(I_forward)
        conv[k+1] = (abs((Refb-Refa)/Refa) )
        if conv[k+1] > 10
            conv[k+1] = 10
        end

        println("convergence: $(conv[k+1])")
        if k > 3
            if conv[k] < tol && conv[k+1] < tol
                break
            end
        end     
        Refa = Refb

    end



    figure()
    plot(λ*1e9, I_forward[:,end])
    plot(λ*1e9, I_backward[:,1])
    yscale("log")

end

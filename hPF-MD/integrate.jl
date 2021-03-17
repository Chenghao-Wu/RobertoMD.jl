

export BerendsenThermostat,calc_totalkineticenergy,calc_inst_temp,apply_1st_integration!,apply_2st_integration!

function calc_totalkineticenergy(atoms::Vector{atom},velocities::Vector{velocity})
    K_tot::Float64=0
    for velocityii=1:length(velocities)
        K_tot+=sum(velocities[velocityii].velocity.^2 .* 0.5 .* atoms[velocityii].mass)
    end
    return K_tot
end

function calc_inst_temp(atoms::Vector{atom},totalkineticenergy::Float64)
    dim::Int64=3
    return totalkineticenergy*2/(length(atoms)*dim-dim)
end

function apply_1st_integration!(args::system,atoms::Vector{atom},velocities::Vector{velocity},forces::Vector{force},::NoThermostat) end

function apply_2nd_integration!(args::system,atoms::Vector{atom},velocities::Vector{velocity},forces::Vector{force},::NoThermostat) end


struct BerendsenThermostat <: Thermostat 
    τ::Float64
end

function apply_1st_integration!(args::system,atoms::Vector{atom},velocities::Vector{velocity},forces::Vector{force},Thermostat::BerendsenThermostat) 
    # update positions
    for posii=1:length(atoms)
        atoms[posii].pos .+= (velocities[posii].velocity .+ 0.5 .* forces[posii].force ./ atoms[posii].mass .* args.dt) .* args.dt
    end 

    # update velocities
    for velii=1:length(atoms)
        velocities[velii].velocity .+= 0.5 .* forces[velii].force ./ atoms[velii].mass .* args.dt
    end 

    # update images
    for posii=1:length(atoms)
        atoms[posii].image .+= RInt.( atoms[posii].pos ./ args.box)
        atoms[posii].pos .-= args.box .* RInt.( atoms[posii].pos ./ args.box)
    end 
end

function apply_2nd_integration!(args::system,atoms::Vector{atom},velocities::Vector{velocity},forces::Vector{force},Thermostat::BerendsenThermostat) 

    totalkineticenergy=calc_totalkineticenergy(atoms,velocities)
    curr_T=calc_inst_temp(atoms,totalkineticenergy)
    scalT  = sqrt(1.0+args.dt*(args.temp/curr_T-1.0)/Thermostat.τ)
    
    # update regulated velocities
    for velii=1:length(atoms)
        velocities[velii].velocity = (velocities[velii].velocity .+0.5 .* forces[velii].force ./ atoms[velii].mass .* args.dt) .* scalT
    end 

end



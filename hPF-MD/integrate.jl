

export BerendsenThermostat,AndersenThermostat

function calc_totalkineticenergy(atoms::Vector{atom},velocities::Vector{velocity})
    K_tot::Float64=0
    for velocityii=1:length(velocities)
        for i=1:3
            K_tot+=sum(velocities[velocityii].velocity[i]^2 * 0.5 * atoms[velocityii].mass)
        end
    end
    return K_tot
end

function calc_inst_temp(atoms::Vector{atom},totalkineticenergy::Float64)
    dim::Int64=3
    T::Float64=totalkineticenergy*2/(length(atoms)*dim-dim)
    return T
end

function apply_1st_integration!(args::system,atoms::Vector{atom},velocities::Vector{velocity},forces::Vector{force},::NoThermostat) end

function apply_2nd_integration!(args::system,atoms::Vector{atom},velocities::Vector{velocity},forces::Vector{force},::NoThermostat) end


struct BerendsenThermostat <: Thermostat 
    τ::Float64
end

function apply_1st_integration!(args::system,atoms::Vector{atom},velocities::Vector{velocity},forces::Vector{force},Thermostat::BerendsenThermostat) 
    # update positions
    for posii=1:length(atoms)
        atoms[posii].pos[1] += (velocities[posii].velocity[1] + 0.5 * forces[posii].force[1] / atoms[posii].mass * args.dt) * args.dt
        atoms[posii].pos[2] += (velocities[posii].velocity[2] + 0.5 * forces[posii].force[2] / atoms[posii].mass * args.dt) * args.dt
        atoms[posii].pos[3] += (velocities[posii].velocity[3] + 0.5 * forces[posii].force[3] / atoms[posii].mass * args.dt) * args.dt
    end 

    # update velocities
    for velii=1:length(atoms)
        velocities[velii].velocity[1] += 0.5 * forces[velii].force[1] / atoms[velii].mass * args.dt
        velocities[velii].velocity[2] += 0.5 * forces[velii].force[2] / atoms[velii].mass * args.dt
        velocities[velii].velocity[3] += 0.5 * forces[velii].force[3] / atoms[velii].mass * args.dt
    end 

    # update images
    for posii=1:length(atoms)
        atoms[posii].image[1] += fld( atoms[posii].pos[1] , args.box[1])
        atoms[posii].image[2] += fld( atoms[posii].pos[2] , args.box[2])
        atoms[posii].image[3] += fld( atoms[posii].pos[3] , args.box[3])
        atoms[posii].pos[1]   -= args.box[1] * fld( atoms[posii].pos[1] , args.box[1])
        atoms[posii].pos[2]   -= args.box[2] * fld( atoms[posii].pos[2] , args.box[2])
        atoms[posii].pos[3]   -= args.box[3] * fld( atoms[posii].pos[3] , args.box[3])
    end 
end

function apply_2nd_integration!(args::system,atoms::Vector{atom},velocities::Vector{velocity},forces::Vector{force},Thermostat::BerendsenThermostat) 

    totalkineticenergy=calc_totalkineticenergy(atoms,velocities)
    curr_T = calc_inst_temp(atoms,totalkineticenergy)
    scalT  = sqrt(1.0+args.dt*(args.temp/curr_T-1.0)/Thermostat.τ)
    
    # update regulated velocities
    for velii=1:length(atoms)
        velocities[velii].velocity[1] = (velocities[velii].velocity[1] +0.5 * forces[velii].force[1] / atoms[velii].mass * args.dt) * scalT
        velocities[velii].velocity[2] = (velocities[velii].velocity[2] +0.5 * forces[velii].force[2] / atoms[velii].mass * args.dt) * scalT
        velocities[velii].velocity[3] = (velocities[velii].velocity[3] +0.5 * forces[velii].force[3] / atoms[velii].mass * args.dt) * scalT
    end 

end

struct AndersenThermostat <: Thermostat 
    σ::Float64
end

function apply_1st_integration!(args::system,atoms::Vector{atom},velocities::Vector{velocity},forces::Vector{force},Thermostat::AndersenThermostat) 
    # update positions
    for posii=1:length(atoms)
        atoms[posii].pos[1] += (velocities[posii].velocity[1] + 0.5 * forces[posii].force[1] / atoms[posii].mass * args.dt) * args.dt
        atoms[posii].pos[2] += (velocities[posii].velocity[2] + 0.5 * forces[posii].force[2] / atoms[posii].mass * args.dt) * args.dt
        atoms[posii].pos[3] += (velocities[posii].velocity[3] + 0.5 * forces[posii].force[3] / atoms[posii].mass * args.dt) * args.dt
    end 

    # update velocities
    for velii=1:length(atoms)
        velocities[velii].velocity[1] += 0.5 * forces[velii].force[1] / atoms[velii].mass * args.dt
        velocities[velii].velocity[2] += 0.5 * forces[velii].force[2] / atoms[velii].mass * args.dt
        velocities[velii].velocity[3] += 0.5 * forces[velii].force[3] / atoms[velii].mass * args.dt
    end 

    # update images
    for posii=1:length(atoms)
        atoms[posii].image[1] += fld( atoms[posii].pos[1] , args.box[1])
        atoms[posii].image[2] += fld( atoms[posii].pos[2] , args.box[2])
        atoms[posii].image[3] += fld( atoms[posii].pos[3] , args.box[3])
        atoms[posii].pos[1]   -= args.box[1] * fld( atoms[posii].pos[1] , args.box[1])
        atoms[posii].pos[2]   -= args.box[2] * fld( atoms[posii].pos[2] , args.box[2])
        atoms[posii].pos[3]   -= args.box[3] * fld( atoms[posii].pos[3] , args.box[3])
    end 
end

function apply_2nd_integration!(args::system,atoms::Vector{atom},velocities::Vector{velocity},forces::Vector{force},Thermostat::AndersenThermostat) 
    # update regulated velocities
    for velii=1:length(atoms)
        d = Normal(0.0, sqrt(args.temp/atoms[velii].mass))
        rn=rand(1)
        if rn[1] < Thermostat.σ * args.dt
            velocities[velii].velocity = rand(d,3)#randn(rng,Float64) * sqrt(args.temp/atoms[velii].mass)
        else
            velocities[velii].velocity[1] += 0.5 * forces[velii].force[1] / atoms[velii].mass * args.dt
            velocities[velii].velocity[2] += 0.5 * forces[velii].force[2] / atoms[velii].mass * args.dt
            velocities[velii].velocity[3] += 0.5 * forces[velii].force[3] / atoms[velii].mass * args.dt
        end
    end
end
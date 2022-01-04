

function Thermostat!(sys::System,comm::MPI.Comm,root::Int64)
    if sys.thermostat!=NoThermostat
        apply_thermostat!(sys,sys.thermostat,comm,root)
    end
end

function apply_thermostat!(sys::System,thermost::BerendsenNVT,comm::MPI.Comm,root::Int64)
    if sys.current_temp[1]==0.0
        sys.current_temp[1]=1e-8
    end
    scalT  = sqrt(1.0+sys.dt*(sys.temp/sys.current_temp[1]-1.0)/thermost.τ)
    apply_thermostat_particle!.(sys.vels,sys.forces,sys.masses,scalT,sys.dt,Ref(thermost))
    #@show sys.forces
end

function apply_thermostat_particle!(vels::Array{Float64,1},
                                    forces::Array{Float64,1},
                                    masses::Float64,
                                    scalT::Float64,
                                    dt::Float64,
                                    thermost::BerendsenNVT)
    vels[1] = vels[1] * scalT
    vels[2] = vels[2] * scalT
    vels[3] = vels[3] * scalT
end

function apply_thermostat!(sys::System,thermost::LangevinNVT,comm::MPI.Comm,root::Int64)
    coeff=sqrt(6.0 * thermost.gamma * sys.temp*thermost.invdt)
    apply_thermostat_particle!.(sys.vels,sys.forces,sys.masses,coeff,sys.dt,Ref(thermost))
end

function apply_thermostat_particle!(vels::Array{Float64,1},
                                    forces::Array{Float64,1},
                                    masses::Float64,
                                    coeff::Float64,
                                    dt::Float64,
                                    thermost::LangevinNVT)
    r_force_x = RandomNumber()*coeff - thermost.gamma*vels[1]
    r_force_y = RandomNumber()*coeff - thermost.gamma*vels[2]
    r_force_z = RandomNumber()*coeff - thermost.gamma*vels[3]
    vels[1] += 0.5 * dt * r_force_x / masses
    vels[2] += 0.5 * dt * r_force_y / masses
    vels[3] += 0.5 * dt * r_force_z / masses
    forces[1] += r_force_x
    forces[2] += r_force_y
    forces[3] += r_force_z
end
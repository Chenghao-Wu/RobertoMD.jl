

function Thermostat!(sys::System,comm::MPI.Comm,root::Int64)
    if sys.thermostat!=NoThermostat
        apply_thermostat!(sys,sys.thermostat,comm,root)
    end
end

function apply_thermostat!(sys::System,thermost::BerendsenNVT,comm::MPI.Comm,root::Int64)
    if sys.current_temp==0.0
        sys.current_temp=1e-8
    end
    scalT  = sqrt(1.0+sys.dt*(sys.temp/sys.current_temp-1.0)/thermost.τ)
    for atomii=1:sys.local_N
        @fastmath @inbounds begin
        sys.vels[atomii,1] = sys.vels[atomii,1] .* scalT
        sys.vels[atomii,2] = sys.vels[atomii,2] .* scalT
        sys.vels[atomii,3] = sys.vels[atomii,3] .* scalT
        end
    end 
end

function apply_thermostat!(sys::System,thermost::LangevinNVT,comm::MPI.Comm,root::Int64)
    coeff=sqrt(6.0 * thermost.gamma * sys.temp*thermost.invdt)
    for atomii=1:sys.local_N
        @fastmath @inbounds begin
        #Random.seed!(MPI.Comm_rank(comm)*sys.local_N+atomii)
        r_force_x = RandomNumber()*coeff - thermost.gamma*sys.vels[atomii,1]
        r_force_y = RandomNumber()*coeff - thermost.gamma*sys.vels[atomii,2]
        r_force_z = RandomNumber()*coeff - thermost.gamma*sys.vels[atomii,3]
        sys.vels[atomii,1] += 0.5 * sys.dt * r_force_x / sys.masses[atomii]
        sys.vels[atomii,2] += 0.5 * sys.dt * r_force_y / sys.masses[atomii]
        sys.vels[atomii,3] += 0.5 * sys.dt * r_force_z / sys.masses[atomii]
        sys.forces[atomii,1] += r_force_x
        sys.forces[atomii,2] += r_force_y
        sys.forces[atomii,3] += r_force_z
        end
    end 
end
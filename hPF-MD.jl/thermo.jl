function UpdateThermo!(sys::System,comm::MPI.Comm,root::Int64)
    UpdateTemperature!(sys,comm,root)
    UpdateMomentum!(sys,comm,root)
end

function kinetic_energy(vel::Array{Float64,1},masses::Float64)
    k=vel[1]*vel[1]*0.5*masses+vel[2]*vel[2]*0.5*masses+vel[3]*vel[3]*0.5*masses
    return k
end

function UpdateTemperature!(sys::System,comm::MPI.Comm,root::Int64)
    Ktot=0.0
    Ktot_local = sum(kinetic_energy.(sys.vels,sys.masses))
    MPI.Barrier(comm)
    Ktot=MPI.Allreduce(Ktot_local, +, comm)
    MPI.Barrier(comm)
    temp=Ktot*2/(sys.total_N*sys.dim-sys.dim)
    sys.current_temp[1]=temp
end

function UpdateMomentum!(sys::System,comm::MPI.Comm,root::Int64)
    momentum_local = sum(sys.vels .* sys.masses)
    momentum_tot=MPI.Allreduce(momentum_local, +, comm)
    momentum_particle=momentum_tot/(sys.total_N*sys.dim-sys.dim)
    sys.current_momentum[1]=momentum_particle[1]+momentum_particle[2]+momentum_particle[3]
end
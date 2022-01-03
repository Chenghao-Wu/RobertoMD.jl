function UpdateThermo!(sys::System,comm::MPI.Comm,root::Int64)
    UpdateTemperature!(sys,comm,root)
    UpdateMomentum!(sys,comm,root)
end


function UpdateTemperature!(sys::System,comm::MPI.Comm,root::Int64)
    Ktot=0.0
    Ktot_local = sum(sys.vels.^2 .* 0.5 .* sys.masses)
    MPI.Barrier(comm)
    Ktot=MPI.Allreduce(Ktot_local, +, comm)
    MPI.Barrier(comm)
    temp=Ktot*2/(sys.total_N*sys.dim-sys.dim)
    sys.current_temp=temp
    if MPI.Comm_rank(comm)==root
        if sys.first_step || sys.current_step%sys.thermofreq==0.0
            @show temp
        end
    end
end

function UpdateMomentum!(sys::System,comm::MPI.Comm,root::Int64)
    momentum_local = sum(sys.vels .* sys.masses)
    momentum_tot=MPI.Allreduce(momentum_local, +, comm)
    if MPI.Comm_rank(comm)==root
        if sys.first_step || sys.current_step%sys.thermofreq==0.0
            @show momentum_tot
        end
    end
end
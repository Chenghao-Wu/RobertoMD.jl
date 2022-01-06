export Simulate


function Initialize_MD!(sys::System,comm::MPI.Comm,root::Int64)
    CalculateForce!(sys,comm)
    UpdateThermo!(sys,comm,root)
    DumpInformation(sys,comm,root)
    LoggerInformation!(sys,comm,root)
    sys.first_step[1]=false
end

function Initialize_LocalSystem!(sys,
                                local_types,
                                local_masses,
                                local_energies,
                                local_coords,
                                local_forces,
                                local_vels,
                                local_bonds)
        for i in 1:size(local_types)[1]
            push!(sys.atoms.types,local_types[i])
        end

        for i in 1:size(local_masses)[1]
            push!(sys.atoms.masses,local_masses[i])
        end

        for i in 1:size(local_energies)[1]
            push!(sys.atoms.energies,local_energies[i])
        end

        for i in 1:size(local_coords)[1]
            push!(sys.atoms.coords,local_coords[i,1:3])
        end

        for i in 1:size(local_forces)[1]
            push!(sys.atoms.forces,local_forces[i,1:3])
        end

        for i in 1:size(local_vels)[1]
            push!(sys.atoms.vels,local_vels[i,1:3])
        end

        for i in 1:size(local_bonds)[1]
            push!(sys.atoms.bonds,local_bonds[i,1:3])
        end
        if sys.field!=NoField()
            for i in 1:size(local_types)[1]
                push!(sys.field.particle2cell,zeros(8))
                push!(sys.field.cellindex,zeros(Int64,8,3))
            end
        end

end


function Simulate(inputs::Dict,config::Dict)
    
    

    MPI.Init()
    comm=MPI.COMM_WORLD
    root=0
    my_rank = MPI.Comm_rank(comm)

    # Data Structure

    if my_rank==root
        sys=System(inputs,config,comm,root)
    else
        sys=System(Dict(),Dict(),comm,root)
    end
    if my_rank==root
        p = Progress(sys.steps; showspeed=true)
    end
    
    sys = MPI.bcast(sys, root, comm)
    
    local_types = MPI.Scatterv(sys.atoms.types_all, sys.mpi.commsizes, root, comm)
    local_masses = MPI.Scatterv(sys.atoms.masses_all, sys.mpi.commsizes, root, comm)
    local_energies = MPI.Scatterv(sys.atoms.energies_all, sys.mpi.commsizes, root, comm)
    
    local_coords = DeSerializeArray(MPI.Scatterv(SerializeArray(sys.atoms.coords_all), sys.mpi.commsizes*3, root, comm))
    local_forces = DeSerializeArray(MPI.Scatterv(SerializeArray(sys.atoms.forces_all), sys.mpi.commsizes*3, root, comm))
    local_vels = DeSerializeArray(MPI.Scatterv(SerializeArray(sys.atoms.vels_all), sys.mpi.commsizes*3, root, comm))
     
    if sys.bdinter!=BdInter[]
        local_bonds = DeSerializeBondArray(MPI.Scatterv(SerializeArray(sys.atoms.bonds_all), sys.mpi.bond_commsizes*3, root, comm))
    else
        local_bonds=zeros(0,3)
    end
    
    Initialize_LocalSystem!(sys,
                            local_types,
                            local_masses,
                            local_energies,
                            local_coords,
                            local_forces,
                            local_vels,
                            local_bonds)

    # Main Algorithm    
    Initialize_MD!(sys,comm,root)

    
    for stepi in sys.current_step[1]+1:sys.current_step[1]+sys.steps-1
        sys.current_step[1]=stepi
        UpdatePosition!(sys,sys.integrator)
        UpdateVelocity!(sys,sys.integrator,comm)
        MPI.Barrier(comm)
        ClearForce!(sys)
        CalculateForce!(sys,comm)
        MPI.Barrier(comm)
        UpdateThermo!(sys,comm,root)
        MPI.Barrier(comm)
        UpdateVelocity!(sys,sys.integrator,comm)
        MPI.Barrier(comm)
        Thermostat!(sys,comm,root)
        MPI.Barrier(comm)
        if my_rank==root
            ProgressMeter.next!(p)
        end
        UpdateThermo!(sys,comm,root)
        LoggerInformation!(sys,comm,root)
        DumpInformation(sys,comm,root)
        RestartInformation(sys,comm,root)
    end
    MPI.Finalize()
end
 

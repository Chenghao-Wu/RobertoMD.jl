export Simulate


function Initialize_MD!(sys::System,comm::MPI.Comm,root::Int64)
    CalculateForce!(sys,comm)
    UpdateThermo!(sys,comm,root)
    #DumpInformation(sys,comm,root)
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
            push!(sys.types,local_types[i])
        end

        for i in 1:size(local_masses)[1]
            push!(sys.masses,local_masses[i])
        end

        for i in 1:size(local_energies)[1]
            push!(sys.energies,local_energies[i])
        end

        for i in 1:size(local_coords)[1]
            push!(sys.coords,local_coords[i,1:3])
        end

        for i in 1:size(local_forces)[1]
            push!(sys.forces,local_forces[i,1:3])
        end

        for i in 1:size(local_vels)[1]
            push!(sys.vels,local_vels[i,1:3])
        end

        for i in 1:size(local_bonds)[1]
            push!(sys.bonds,local_bonds[i,1:3])
        end
        if sys.field!=NoField()
            for i in 1:size(local_types)[1]
                push!(sys.particle2cell,zeros(8))
                push!(sys.cellindex,zeros(Int64,8,3))
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
        sys=System(inputs,config,comm)
    else
        sys=System(Dict(),Dict(),comm)
    end
    if my_rank==root
        p = Progress(sys.steps; showspeed=true)
    end
    
    sys = MPI.bcast(sys, root, comm)
    
    local_types = MPI.Scatterv(sys.types_all, sys.commsizes, root, comm)
    local_masses = MPI.Scatterv(sys.masses_all, sys.commsizes, root, comm)
    local_energies = MPI.Scatterv(sys.energies_all, sys.commsizes, root, comm)
    
    local_coords = DeSerializeArray(MPI.Scatterv(SerializeArray(sys.coords_all), sys.commsizes*3, root, comm))
    local_forces = DeSerializeArray(MPI.Scatterv(SerializeArray(sys.forces_all), sys.commsizes*3, root, comm))
    local_vels = DeSerializeArray(MPI.Scatterv(SerializeArray(sys.vels_all), sys.commsizes*3, root, comm))
    
    if sys.bdinter!=NoBond()
        local_bonds = DeSerializeBondArray(MPI.Scatterv(SerializeArray(sys.bonds_all), sys.bond_commsizes*3, root, comm))
    else
        local_bonds=zeros(Int64,0,0)
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
    #
    for stepi in sys.current_step[1]+1:sys.current_step[1]+sys.steps
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
        #DumpInformation(sys,comm,root)
        if my_rank==root
            ProgressMeter.next!(p)
        end
    end
    MPI.Finalize()
end
 

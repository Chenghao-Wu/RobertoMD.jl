export Simulate


function Initialize!(sys::System,comm::MPI.Comm,root::Int64)
    CalculateForce!(sys,comm)
    #UpdateThermo!(sys,comm,root)
    sys.first_step=false
end

function Simulate(inputs::Dict,config::Dict)

    MPI.Init()
    comm=MPI.COMM_WORLD
    root=0
    my_rank = MPI.Comm_rank(comm)

    # Data Structure

    if my_rank==root
        sys=System(inputs,config,comm)
        configuration=Configuration(config)
        balancingMPI=BalancingMPI(configuration,MPI.Comm_size(comm))
        # Initialize Internal Data Structure
        
        sys.pbc=PBC(;xhi=configuration.box[1],yhi=configuration.box[2],zhi=configuration.box[3])
        sys.rho0=sum(balancingMPI.commsizes)/(sys.pbc.Lx^3)
        sys.commsizes=balancingMPI.commsizes
        sys.bond_commsizes=balancingMPI.bond_commsizes
        
        types=InitTypes(configuration,balancingMPI)
        sys.num_atomtype=length(unique(types))
        sys.total_N=length(types)
        masses=InitMasses(configuration,balancingMPI)
        coordinates=InitCoords(configuration,balancingMPI)
        forces=InitForces(configuration)
        energies=InitEnergies(configuration)
        vels=InitVels(configuration,sys.velocity)
        if sys.bdinter!=NoBond
            bonds=InitBonds(configuration,balancingMPI)
        end
        if sys.field!=NoField()

            if "uniform mesh" in keys(inputs["Canonical field"])
                sys.mesh=UniformMesh(sys.pbc,inputs["Canonical field"]["Lcell"])
            end

            sys.local_field=[zeros(sys.mesh.num_cells*sys.mesh.num_cells*sys.mesh.num_cells) for i in 1:sys.num_atomtype]
            sys.global_field=[zeros(sys.mesh.num_cells*sys.mesh.num_cells*sys.mesh.num_cells) for i in 1:sys.num_atomtype]
            sys.field_gradient=[zeros(sys.mesh.num_cells*sys.mesh.num_cells*sys.mesh.num_cells*3) for i in 1:sys.num_atomtype]
        end
    else
        sys=System(Dict(),Dict(),comm)
        forces=zeros(0,0)
        coordinates=zeros(0,0)
        energies=zeros(0)
        types=zeros(Int64,0)
        masses=zeros(Float64,0)
        bonds=zeros(Int64,0,0)
        vels=zeros(0,0)
    end

    sys = MPI.bcast(sys, root, comm)
    
    sys.types = MPI.Scatterv(types, sys.commsizes, root, comm)
    sys.masses = MPI.Scatterv(masses, sys.commsizes, root, comm)
    sys.energies = MPI.Scatterv(energies, sys.commsizes, root, comm)
    
    sys.coords = DeSerializeArray(MPI.Scatterv(SerializeArray(coordinates), sys.commsizes*3, root, comm))
    sys.forces = DeSerializeArray(MPI.Scatterv(SerializeArray(forces), sys.commsizes*3, root, comm))
    sys.vels = DeSerializeArray(MPI.Scatterv(SerializeArray(vels), sys.commsizes*3, root, comm))
    
    if sys.bdinter!=NoBond()
        sys.bonds = DeSerializeBondArray(MPI.Scatterv(SerializeArray(bonds), sys.bond_commsizes*3, root, comm))
    end
    sys.local_N=size(sys.types)[1]
    if sys.field!=NoField()
        sys.particle2cell=zeros(sys.local_N,8)
        sys.cellindex=[zeros(Int64,8,3) for i in 1:sys.local_N]
    end
    
    # Main Algorithm
    Initialize!(sys,comm,root)

    for stepi in sys.current_step+1:sys.current_step+sys.steps
        sys.current_step=stepi
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
    end
    MPI.Finalize()
end
 

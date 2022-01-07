export System,ReadInput

function ReadInput(filename::String)
    input = open(filename, "r")
    inputs=JSON.parse(input)
    close(input)
    return inputs
end

struct System
    # User Inputs Simulation Setup
    dt::Float64
    temp::Float64 
    steps::Int64
    start_step::Int64
    atoms::Atoms
    thermostat::Array{Thermostat,1}
    integrator::Integrator
    bdinter::Array{BdInter,1}
    field::Field
    velocity::Velocity
    thermologger::Logger
    trjdump::Dump
    restart::Restart
    pbc::BC

    first_step::Array{Bool,1}
    current_step::Array{Int64,1}
    current_temp::Array{Float64,1}
    current_momentum::Array{Float64,1}
    mpi::BalancingMPI
    
    # Initialize Constructor
    function System(input::Dict,config::Dict,comm::MPI.Comm,root::Int64)

        # local worker ID 
        my_rank = MPI.Comm_rank(comm)
        comm_size=MPI.Comm_size(comm)

        # Initialize Simulation Setup
        start_step=1
        dt=0.005
        temp=1.0
        steps=100

        atoms=Atoms()
        pbc=NoPBC()

        bdinter=Array{BdInter,1}(undef,0)
        field=NoField()
        integrator=NoIntegrator()
        thermostat=Array{Thermostat,1}(undef,0)
        velocity=NoVelocity()
        trjdump=NoDump()
        restart=NoRestart()

        first_step=[true]
        current_step=[1]
        current_temp=[0.0]
        current_momentum=[0.0]

        thermologger=NoLogger()

        configuration=Configuration(Dict())
        mpi=BalancingMPI(configuration,comm_size)
        
        if my_rank==root

            if "dt" in keys(input)
                dt=input["dt"]
                @info "Timestep dt: $(dt)"
            end
            
            if "temp" in keys(input)
                temp=input["temp"]
                @info "Temperature T: $(temp)"
            end
            if "steps" in keys(input)
                steps=input["steps"]
                @info "Simulation Total Steps: $(steps)"
            end
            configuration=Configuration(config)
            mpi=BalancingMPI(configuration,MPI.Comm_size(comm))

            start_step=configuration.start
            current_step[1]=configuration.start

            pbc=PBC(;xhi=configuration.box[1],yhi=configuration.box[2],zhi=configuration.box[3])
            atoms.rho0[1]=sum(mpi.commsizes)/(pbc.Lx^3)
            
            types_all=InitTypes(configuration,mpi)

            atoms.num_atomtype[1]=length(unique(types_all))
            atoms.total_N[1]=length(types_all)

            if "velocity verlet" in keys(input)
                integrator=VVIntegrator(dt)
                @info "Velocity Velert Integrator"
            end

            if "gaussian velocity" in keys(input)
                temp=input["gaussian velocity"]["T"]
                seed=input["gaussian velocity"]["Seed"]
                velocity=GaussianVelocity(temp,seed)
                @info "Gaussian Velocity Initialization with Temperature: $(temp) and Seed:  $(seed)"
            elseif "zero velocity" in keys(input)
                velocity=ZeroVelocity()
                @info "Zero Velocity Initialization"
            end
            if "restart" in keys(input)
                if input["restart"]["JSONRestart"]
                    filename=input["restart"]["file"]
                    freq=input["restart"]["freq"]
                    if freq==-1
                        freq=Int(1e10)
                    end
                    restart=JSONRestart(filename,freq)
                    
                end
            end

            masses_all=InitMasses(configuration,mpi)
            coords_all=InitCoords(configuration,mpi)
            forces_all=InitForces(configuration)
            vels_all=InitVels(configuration,mpi,velocity)
            energies_all=InitEnergies(configuration)

            init_atoms!(atoms,
                        types_all,
                        masses_all,
                        coords_all,
                        forces_all,
                        vels_all,
                        energies_all)

            if "harmonic bond" in keys(input)
                bonds_all=InitBonds(configuration,mpi)
                init_bonds!(atoms,bonds_all)
            end

            if "thermo information" in keys(input)
                thermologger=init(NoLogger(),input,start_step)
            end

            if "BerendsenNVT" in keys(input)
                thermostat_=init(BerendsenNVT(),input)
                push!(thermostat,thermostat_)
                @info "BerendsenNVT Thermostat with tau $(thermostat_.τ)"
            elseif "LangevinNVT" in keys(input)
                thermostat_=init(LangevinNVT(),input)
                push!(thermostat,thermostat_)
                @info "LangevinNVT Thermostat with gamma $(thermostat_.gamma)"
            end

            if "harmonic bond" in keys(input)
                bdinter_=init(Harmonic(),input)
                push!(bdinter,bdinter_)
                @info "Employing Harmonic Bond Potential with K= $(input["harmonic bond"]["k"]) and l0= $(input["harmonic bond"]["l0"])"
            end

            if "Canonical field" in keys(input)
                field=init(CanonicalField(),input,configuration)
                @info "Employing Canonical Field with Cell Length $(input["Canonical field"]["Lcell"]), Field Update every $(input["Canonical field"]["update"])"
            end

            if "LAMMPSTrj" in keys(input)
                filename=input["LAMMPSTrj"]["file"]
                freq=input["LAMMPSTrj"]["freq"]
                trjdump=LAMMPSTrj(filename,freq)
            end
        end
        

        new(# User Inputs Simulation Setup
        dt::Float64,
        temp::Float64 ,
        steps::Int64,
        start_step::Int64,
    
        atoms::Atoms,
        thermostat::Array{Thermostat,1},
        integrator::Integrator,
        bdinter::Array{BdInter,1},
        field::Field,
        velocity::Velocity,
        thermologger::Logger,
        trjdump::Dump,
        restart::Restart,
        pbc::BC,
    
        first_step::Array{Bool,1},
        current_step::Array{Int64,1},
        current_temp::Array{Float64,1},
        current_momentum::Array{Float64,1},
        mpi::BalancingMPI
        )
    end
end

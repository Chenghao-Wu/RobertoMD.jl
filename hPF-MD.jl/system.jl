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
    thermofreq::Int64
    thermostat::Thermostat
    integrator::Integrator
    bdinter::BdInter
    field::Field
    mesh::Mesh
    velocity::Velocity
    first_step::Array{Bool,1}
    current_step::Array{Int64,1}
    current_temp::Array{Float64,1}
    dim::Int64
    trjdump::Dump
    commsizes::Array{Int64,1}
    bond_commsizes::Array{Int64,1}

    num_atomtype::Int64
    pbc::BC
    rho0::Float64
    local_field::Array{Array{Float64,1},1}
    global_field::Array{Array{Float64,1},1}
    field_gradient::Array{Array{Float64,1},1}
    particle2cell::Array{Array{Float64,1},1}
    cellindex::Array{Array{Int64,2},1}

    total_N::Int64
    local_N::Int64
    
    types::Array{Int64,1}
    masses::Array{Float64,1}
    coords::Vector{Array{Float64,1}}
    vels::Vector{Array{Float64,1}}
    forces::Vector{Array{Float64,1}}
    energies::Array{Float64,1}
    bonds::Vector{Array{Int64,1}}

    normal::Normal{Float64}

    coords_all::Array{Float64,2}
    types_all::Array{Int64,1}
    masses_all::Array{Float64,1}
    vels_all::Array{Float64,2}
    forces_all::Array{Float64,2}
    energies_all::Array{Float64,1}
    bonds_all::Array{Int64,2}

    # Initialize Constructor
    function System(input::Dict,config::Dict,comm::MPI.Comm)
        # define the root worker
        root = 0
        # local worker ID 
        my_rank = MPI.Comm_rank(comm)
        comm_size=MPI.Comm_size(comm)

        # Initialize Simulation Setup
        dt=0.005
        temp=1.0
        steps=100
        thermofreq=1
        bdinter=NoBond()
        field=NoField()
        integrator=NoIntegrator()
        thermostat=NoThermostat()
        velocity=NoVelocity()
        trjdump=NoDump()

        commsizes=zeros(Int64,0)
        bond_commsizes=zeros(Int64,0)

        dim=3
        total_N=0
        num_atomtype=0
        local_N=0
        first_step=[true]
        current_step=[1]
        current_temp=[0.0]

        local_field=[zeros(0)]
        global_field=[zeros(0)]
        field_gradient=[zeros(0)]
        
        particle2cell=Array{Array{Float64,1},1}(undef,0)
        cellindex=Array{Array{Int64,2},1}(undef,0)

        normal=Normal()

        rho0=0.0

        mesh=NoMesh()
        pbc=NoPBC()

        configuration=Configuration(Dict())
        balancingMPI=BalancingMPI(configuration,comm_size)
        
        types=Vector{Int64}(undef,0)
        masses=Vector{Float64}(undef,0)
        coords=Array{Array{Float64,1},1}(undef,0)
        vels=Array{Array{Float64,1},1}(undef,0)
        forces=Array{Array{Float64,1},1}(undef,0)
        energies=Vector{Float64}(undef,0)
        bonds=Vector{Vector{Int64}}(undef,0)
        
        coords_all=zeros(Float64,0,0)
        types_all=zeros(Int64,0)
        masses_all=zeros(Float64,0)
        vels_all=zeros(Float64,0,0)
        forces_all=zeros(Float64,0,0)
        energies_all=zeros(Float64,0)
        bonds_all=zeros(Int64,0,0)

        if my_rank==root
            
            if "dt" in keys(input)
                dt=input["dt"]
                @info "Timestep dt: $(dt)"
            end
            if "start" in keys(input)
                current_step[1]=input["start"]
                @info "Start Step dt: $(current_step)"
            end
            if "temp" in keys(input)
                temp=input["temp"]
                @info "Temperature T: $(temp)"
            end
            if "steps" in keys(input)
                steps=input["steps"]
                @info "Simulation Total Steps: $(steps)"
            end
            if "thermofreq" in keys(input)
                thermofreq=input["thermofreq"]
                @info "Log Thermo Information Every: $(thermofreq)"
            end

            if "BerendsenNVT" in keys(input)
                tau=input["BerendsenNVT"]["tau"]
                thermostat=BerendsenNVT(tau)
                @info "BerendsenNVT Thermostat with tau $(tau)"
            elseif "LangevinNVT" in keys(input)
                gamma=input["LangevinNVT"]["gamma"]
                thermostat=LangevinNVT(gamma,dt)
                @info "LangevinNVT Thermostat with gamma $(gamma)"
            end

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

            if "harmonic bond" in keys(input)
                k_matrix=zeros(length(keys(input["harmonic bond"]["k"])))
                l0_matrix=zeros(length(keys(input["harmonic bond"]["l0"])))

                for atomi in 1:length(keys(input["harmonic bond"]["k"]))
                    k_matrix[atomi]=input["harmonic bond"]["k"][string(atomi)]
                    l0_matrix[atomi]=input["harmonic bond"]["l0"][string(atomi)]
                end
                
                bdinter=Harmonic(k_matrix,l0_matrix)
                @info "Employing Harmonic Bond Potential with K= $(input["harmonic bond"]["k"]) and l0= $(input["harmonic bond"]["l0"])"
            end

            if "Canonical field" in keys(input)
                num_atomtypes=length(keys(input["Canonical field"]["χ"]))
                chi=zeros(num_atomtypes,num_atomtypes)
                for atomi in 1:num_atomtypes
                    chi_i=convert(Array{Float64,1},input["Canonical field"]["χ"][string(atomi)])
                    for atomj in 1:num_atomtypes
                        chi[atomi,atomj]=chi_i[atomj]
                    end
                end
                field=CanonicalField(chi,input["Canonical field"]["κ"],input["Canonical field"]["update"],input["Canonical field"]["Lcell"])
                
                @info "Employing Canonical Field with Cell Length $(input["Canonical field"]["Lcell"]), Field Update every $(input["Canonical field"]["update"])"
            end

            if "LAMMPSTrj" in keys(input)
                filename=input["LAMMPSTrj"]["file"]
                freq=input["LAMMPSTrj"]["freq"]
                trjdump=LAMMPSTrj(filename,freq)
            end

            configuration=Configuration(config)
            balancingMPI=BalancingMPI(configuration,MPI.Comm_size(comm))

            pbc=PBC(;xhi=configuration.box[1],yhi=configuration.box[2],zhi=configuration.box[3])
            rho0=sum(balancingMPI.commsizes)/(pbc.Lx^3)
            commsizes=balancingMPI.commsizes
            bond_commsizes=balancingMPI.bond_commsizes
            
            types_all=InitTypes(configuration,balancingMPI)
            num_atomtype=length(unique(types_all))

            total_N=length(types_all)
            masses_all=InitMasses(configuration,balancingMPI)
            coords_all=InitCoords(configuration,balancingMPI)
            forces_all=InitForces(configuration)
            energies_all=InitEnergies(configuration)
            vels_all=InitVels(configuration,velocity)
            if bdinter!=NoBond
                bonds_all=InitBonds(configuration,balancingMPI)
            end
            if field!=NoField()
                if "uniform mesh" in keys(input["Canonical field"])
                    mesh=UniformMesh(pbc,input["Canonical field"]["Lcell"])
                end

                local_field=[zeros(mesh.num_cells*mesh.num_cells*mesh.num_cells) for i in 1:num_atomtype]
                global_field=[zeros(mesh.num_cells*mesh.num_cells*mesh.num_cells) for i in 1:num_atomtype]
                field_gradient=[zeros(mesh.num_cells*mesh.num_cells*mesh.num_cells*3) for i in 1:num_atomtype]

                #particle2cell=zeros(maximum(commsizes),8)
                #cellindex=[zeros(Int64,8,3) for i in 1:maximum(commsizes)]
            end
        end
        

        new(dt::Float64,
            temp::Float64 ,
            steps::Int64,
            thermofreq::Int64,
            thermostat::Thermostat,
            integrator::Integrator,
            bdinter::BdInter,
            field::Field,
            mesh::Mesh,
            velocity::Velocity,
            first_step::Array{Bool,1},
            current_step::Array{Int64,1},
            current_temp::Array{Float64,1},
            dim::Int64,
            trjdump::Dump,
            commsizes::Array{Int64,1},
            bond_commsizes::Array{Int64,1},
        
            num_atomtype::Int64,
            pbc::BC,
            rho0::Float64,
            local_field::Array{Array{Float64,1},1},
            global_field::Array{Array{Float64,1},1},
            field_gradient::Array{Array{Float64,1},1},
            particle2cell::Array{Array{Float64,1},1},
            cellindex::Array{Array{Int64,2},1},
        
            total_N::Int64,
            local_N::Int64,

            types::Array{Int64,1},
            masses::Array{Float64,1},
            coords::Vector{Array{Float64,1}},
            vels::Vector{Array{Float64,1}},
            forces::Vector{Array{Float64,1}},
            energies::Array{Float64,1},
            bonds::Vector{Array{Int64,1}},

            normal::Normal{Float64},

            coords_all::Array{Float64,2},
            types_all::Array{Int64,1},
            masses_all::Array{Float64,1},
            vels_all::Array{Float64,2},
            forces_all::Array{Float64,2},
            energies_all::Array{Float64,1},
            bonds_all::Array{Int64,2})
    end
end

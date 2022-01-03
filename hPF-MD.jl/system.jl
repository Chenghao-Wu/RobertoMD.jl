export System,ReadInput

function ReadInput(filename::String)
    input = open(filename, "r")
    inputs=JSON.parse(input)
    close(input)
    return inputs
end

mutable struct System
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
    first_step::Bool
    current_step::Int64
    current_temp::Float64
    dim::Int64

    commsizes::Array{Int64,1}
    bond_commsizes::Array{Int64,1}

    num_atomtype::Int64
    pbc::BC
    rho0::Float64
    local_field::Array{Array{Float64,1},1}
    global_field::Array{Array{Float64,1},1}
    field_gradient::Array{Array{Float64,1},1}
    particle2cell::Array{Float64,2}
    cellindex::Array{Array{Int64,2},1}

    total_N::Int64
    local_N::Int64
    types::Array{Int64,1}
    masses::Array{Float64,1}
    coords::Array{Float64,2}
    vels::Array{Float64,2}
    forces::Array{Float64,2}
    energies::Array{Float64,1}
    bonds::Array{Int64,2}
    normal::Normal{Float64}
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

        commsizes=zeros(Int64,0)
        bond_commsizes=zeros(Int64,0)

        dim=3
        total_N=0
        num_atomtype=0
        local_N=0
        first_step=true
        current_step=1
        current_temp=0.0

        local_field=[zeros(0)]
        global_field=[zeros(0)]
        field_gradient=[zeros(0)]
        particle2cell=zeros(0,0)
        cellindex=[zeros(Int64,0,0)]

        normal=Normal()

        rho0=0.0

        mesh=NoMesh()
        pbc=NoPBC()

        configuration=Configuration(Dict())
        balancingMPI=BalancingMPI(configuration,comm_size)
        
        if my_rank==root
            
            if "dt" in keys(input)
                dt=input["dt"]
                @info "Timestep dt: $(dt)"
            end
            if "start" in keys(input)
                current_step=input["start"]
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
                bdinter=Harmonic(input["harmonic bond"]["k"],input["harmonic bond"]["l0"])
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
        end
        
        types=zeros(Int64,0)
        masses=zeros(Float64,0)
        coords=zeros(Float64,0,0)
        vels=zeros(Float64,0,0)
        forces=zeros(Float64,0,0)
        energies=zeros(Float64,0)
        bonds=zeros(Int64,0,0)

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
            first_step::Bool,
            current_step::Int64,
            current_temp::Float64,
            dim::Int64,
        
            commsizes::Array{Int64,1},
            bond_commsizes::Array{Int64,1},
        
            num_atomtype::Int64,
            pbc::BC,
            rho0::Float64,
            local_field::Array{Array{Float64,1},1},
            global_field::Array{Array{Float64,1},1},
            field_gradient::Array{Array{Float64,1},1},
            particle2cell::Array{Float64,2},
            cellindex::Array{Array{Int64,2},1},
        
            total_N::Int64,
            local_N::Int64,
            types::Array{Int64,1},
            masses::Array{Float64,1},
            coords::Array{Float64,2},
            vels::Array{Float64,2},
            forces::Array{Float64,2},
            energies::Array{Float64,1},
            bonds::Array{Int64,2},
            normal::Normal{Float64})

    end
end

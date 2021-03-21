export atom,velocity,force,sys_init,run_hPFMD

abstract type Thermostat end

abstract type BondInteraction end
abstract type NonBondInteraction end

abstract type Logger end

abstract type TrajectoryDump end

mutable struct atom
    type::Int64
    mass::Float64
    charge::Float64
    pos::Vector{Float64}
    image::Vector{Int64}
end

mutable struct velocity
    velocity::Vector{Float64}
end

mutable struct force
    force::Vector{Float64}
end

struct bond
    atom_1::Int64
    atom_2::Int64
end

struct angle
    atom_1::Int64
    atom_2::Int64
    atom_3::Int64
end

struct NoThermostat <: Thermostat end
struct NoBondInteractions <: BondInteraction end
struct NoNonBondInteractions <: NonBondInteraction end
struct NoLogger <: Logger end
struct NoTrajectoryDump <: TrajectoryDump end

mutable struct system
    dt::Float64
    temp::Float64 
    steps::Int64
    logfreq::Int64
    logfile::AbstractString
    logmode::AbstractString
    trjdumpfreq::Int64
    trjdumpfile::AbstractString
    trjdumpmode::AbstractString
    box::Vector{Float64}
    data_file::AbstractString
    data_format::AbstractString
    thermostat::Thermostat
    bondinteraction::BondInteraction
    nonbondinteraction::NonBondInteraction
    logger::Logger
    trajectorydump::TrajectoryDump
    firststep::Bool
end

function sys_init()
    dt=0.0
    temp=0.0
    logfreq=0
    steps=0
    box=[0.0,0.0,0.0]
    logfile=""
    logmode="w"
    trjdumpfreq=0
    trjdumpfile=""
    trjdumpmode="w"
    data_file=""
    data_format="LAMMPS Data"
    thermostat=NoThermostat()
    bondinteraction=NoBondInteractions()
    nonbondinteraction=NoNonBondInteractions()
    logger=NoLogger()
    trajectorydump=NoTrajectoryDump()
    firststep=true
    return system(dt,temp,steps,logfreq,logfile,logmode,trjdumpfreq,trjdumpfile,trjdumpmode,box,data_file,data_format,thermostat,bondinteraction,nonbondinteraction,logger,trajectorydump,firststep)
end

function read_data(filename::AbstractString,format::AbstractString)
    #a=Atom(1,1,[1,1,1],[1,1,1])
    universe=Trajectory(filename,'r',format)
    frame = read(universe)
    pos = positions(frame)
    center_pos = [minimum(pos[1,:]),minimum(pos[2,:]),minimum(pos[3,:])] ./ 10 # nanometer
    topo = Topology(frame)
    bonds_in=bonds(topo)
    angles_in=angles(topo)

    number_atoms=length(frame)
    atoms_=Vector{atom}(undef,length(frame))
    number_bonds::Int64=length(bonds_in)/2
    bonds_=Vector{bond}(undef,number_bonds)
    number_angles::Int64=length(angles_in)/3
    angles_=Vector{angle}(undef,number_angles)
    
    img_zero=zeros(Int64,3)
    for atomii =1:number_atoms
        type_=1
        atomi_=Atom(frame,atomii-1) # question about the atom index?
        massi_=mass(atomi_)
        chargei_=charge(atomi_)
        current_pos=pos[1:3,atomii]./10 .- center_pos # nanometer
        atomi=atom(type_,massi_,chargei_,current_pos,img_zero) 
        atoms_[atomii]=atomi
    end

    for bondii =1:number_bonds
        bondi_1=bonds_in[1,bondii] 
        bondi_2=bonds_in[2,bondii] 
        bondi=bond(bondi_1,bondi_2)
        bonds_[bondii]=bondi
    end

    for angleii =1:number_angles
        anglei_1=angles_in[1,angleii] 
        anglei_2=angles_in[2,angleii] 
        anglei_3=angles_in[3,angleii] 
        anglei=angle(anglei_1,anglei_2,anglei_3)
        angles_[angleii]=anglei
    end

    
    velocities_=Vector{velocity}(undef,number_atoms)
    velocity_zero=zeros(3)
    for velocityii=1:number_atoms
        velocity_=velocity(velocity_zero)
        velocities_[velocityii]=velocity_
    end

    force_zero=zeros(3)
    forces_=Vector{force}(undef,number_atoms)
    for forceii=1:number_atoms
        force_=force(force_zero)
        forces_[forceii]=force_
    end

    energy_=zeros(Float64,number_atoms)

    return atoms_,bonds_,angles_,velocities_,forces_,energy_
end

function run_hPFMD(args::system)
    # using the same random seed in all simulations
    Random.seed!(0)

    # read data
    atoms_,bonds_,angles_,velocities_,forces_,energy_=read_data(args.data_file,args.data_format)

    firststep=args.firststep

    # initialize veclocities!
    velociy_init=true
    if velociy_init==false
        velocities_=velocities_
    else
        for velocityii=1:length(velocities_)
            d = Normal(0.0, sqrt(args.temp/atoms_[velocityii].mass))
            velocities_[velocityii]=velocity(rand(d,3))
        end
    end

    logger_file = open(args.logfile, args.logmode)
    trjdump_file = open(args.trjdumpfile, args.trjdumpmode)

    if firststep==true
        clear_forces!(forces_)
        clear_energy!(energy_)

        apply_bonds!(args,atoms_,forces_,energy_,bonds_,args.bondinteraction)
        apply_nonbonds!(args,atoms_,forces_,energy_,args.nonbondinteraction)

        apply_logger!(logger_file,firststep,1,velocities_,energy_,atoms_,args.logger)
        apply_dump(trjdump_file,1,args,atoms_,forces_,velocities_,args.trajectorydump)
        firststep=false
    end

    @showprogress for step_i = 2:args.steps
        
        apply_1st_integration!(args,atoms_,velocities_,forces_,args.thermostat)

        clear_forces!(forces_)
        clear_energy!(energy_)

        apply_bonds!(args,atoms_,forces_,energy_,bonds_,args.bondinteraction)
        apply_nonbonds!(args,atoms_,forces_,energy_,args.nonbondinteraction)

        if step_i%args.logfreq==0.0
            apply_logger!(logger_file,firststep,step_i,velocities_,energy_,atoms_,args.logger)
        end

        if step_i%args.trjdumpfreq==0.0
            apply_dump(trjdump_file,step_i,args,atoms_,forces_,velocities_,args.trajectorydump)
        end

        apply_2nd_integration!(args,atoms_,velocities_,forces_,args.thermostat)
        
    end

end


function InitCoords(configuration::Configuration,balancingMPI::BalancingMPI)
    n_total=sum(configuration.moles_length)
    coords=zeros(n_total,3)
    index=1
    for comm_mole_i in balancingMPI.mole_index_comms
        for mole_i in comm_mole_i
            molecule=configuration.moles[mole_i]
            mole_length=size(molecule.coords)[1]
            coords[index:mole_length+index-1,1:3]=molecule.coords
            index=index+mole_length
        end
    end
    @info "Initialize Coordinates Size $(n_total)"
    return coords
end

function InitMasses(configuration::Configuration,balancingMPI::BalancingMPI)
    n_total=sum(configuration.moles_length)
    masses=zeros(n_total)
    index=1
    for comm_mole_i in balancingMPI.mole_index_comms
        for mole_i in comm_mole_i
            molecule=configuration.moles[mole_i]
            mole_length=size(molecule.coords)[1]
            masses[index:mole_length+index-1]=molecule.mass
            index=index+mole_length
        end
    end
    @info "Initialize Masses Size $(n_total)"
    return masses
end

function InitTypes(configuration::Configuration,balancingMPI::BalancingMPI)
    n_total=sum(configuration.moles_length)
    types=zeros(Int64,n_total)
    index=1
    for comm_mole_i in balancingMPI.mole_index_comms
        for mole_i in comm_mole_i
            molecule=configuration.moles[mole_i]
            mole_length=size(molecule.coords)[1]
            types[index:mole_length+index-1]=molecule.type
            index=index+mole_length
        end
    end
    @info "Initialize Types Size $(n_total)"
    return types
end

function InitBonds(configuration::Configuration,balancingMPI::BalancingMPI)
    n_total=0
    for molei in configuration.moles
        n_total+=size(molei.bonds)[1]
    end
    bonds=zeros(Int64,n_total,3)
    index=1
    for comm_mole_i in balancingMPI.mole_index_comms
        total_atom_index=0
        for mole_i in comm_mole_i
            molecule=configuration.moles[mole_i]
            bond_length=size(molecule.bonds)[1]
            mole_length=size(molecule.coords)[1]
            if bond_length!=0
                bonds[index:bond_length+index-1,1]=molecule.bonds[:,1]
                bonds[index:bond_length+index-1,2:3]=molecule.bonds[:,2:3].+total_atom_index
            end
            total_atom_index+=mole_length
            index=index+bond_length
        end
    end
    @info "Initialize Bonds Size $(n_total)"
    return bonds
end

function InitVels(configuration::Configuration,balancingMPI::BalancingMPI,velocity::ZeroVelocity)
    n_total=sum(configuration.moles_length)
    vels=zeros(n_total,3)
    @info "Initialize Velocities Size $(n_total)"
    return vels
end

function InitVels(configuration::Configuration,balancingMPI::BalancingMPI,velocity::NoVelocity)
    n_total=sum(configuration.moles_length)
    vels=zeros(n_total,3)
    index=1
    for comm_mole_i in balancingMPI.mole_index_comms
        for mole_i in comm_mole_i
            molecule=configuration.moles[mole_i]
            mole_length=size(molecule.vels)[1]
            vels[index:mole_length+index-1,1:3]=molecule.vels
            index=index+mole_length
        end
    end
    @info "Initialize Velocities Size $(n_total)"
    return vels
end

function InitForces(configuration::Configuration)
    n_total=sum(configuration.moles_length)
    forces=zeros(n_total,3)
    @info "Initialize Forces Size $(n_total)"
    return forces
end

function InitEnergies(configuration::Configuration)
    n_total=sum(configuration.moles_length)
    energies=zeros(n_total)
    @info "Initialize Energy Size $(n_total)"
    return energies
end



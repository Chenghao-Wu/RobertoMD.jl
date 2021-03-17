
export TableBondInteractions,OriginalhPFInteractions


function clear_forces!(forces::Vector{force})
    for forceii=1:length(forces)
        forces[forceii].force=[0.0,0.0,0.0]
    end
end

function clear_energy(energy::Vector{Float64})
    return zeros(Float64,length(energy))
end

function apply_bonds!(args::system,atoms::Vector{atom},forces::Vector{force},energy::Vector{Float64},bonds::Vector{bond},::NoBondInteractions) end

struct TableBondInteractions <: BondInteraction 
    potentials::Array{Float64}
    δ::Float64
end

function apply_bonds!(args::system,atoms::Vector{atom},forces::Vector{force},energy::Vector{Float64},bonds::Vector{bond},bondtable::TableBondInteractions) 
    for bondii=1:length(bonds)
        delta_d = atoms[bonds[bondii].atom_1+1].pos .- atoms[bonds[bondii].atom_2+1].pos
        # consider pbc
        delta_d .-= args.box .* RInt.( delta_d ./ args.box)

        r=norm(delta_d)
        index_bond=RInt(r/bondtable.δ)
        if index_bond>200
            print("ERROR!! BOND break")
        end
        force_=bondtable.potentials[index_bond,4]
        energy_=bondtable.potentials[index_bond,3]

        energy[bonds[bondii].atom_1+1] += energy_ * 0.5
        energy[bonds[bondii].atom_2+1] += energy_ * 0.5

        force_divr = force_/r
        
        # calculate forces
        forces[bonds[bondii].atom_1+1].force .+= delta_d .* force_divr
        forces[bonds[bondii].atom_2+1].force .+= -delta_d .* force_divr
    end
end

struct HarmonicBondInteractions <: BondInteraction 
    κ::Float64
    r::Float64
end

function apply_bonds!(args::system,atoms::Vector{atom},forces::Vector{force},energy::Vector{Float64},bonds::Vector{bond},harmonicbond::HarmonicBondInteractions) 
    for bondii=1:length(bonds)
        delta_d = atoms[bonds[bondii].atom_1+1].pos .- atoms[bonds[bondii].atom_2+1].pos

        # consider pbc
        delta_d .-= args.box .* RInt.( delta_d ./ args.box)

        r=norm(delta_d)
        index_bond=RInt(r/bondtable.δ_d)
        force_=bondtable.potentials[index_bond,4]
        force_divr = force_/r
        
        # calculate forces
        forces[bonds[bondii].atom_1+1].force .+= delta_d .* force_divr
        forces[bonds[bondii].atom_2+1].force .+= -delta_d .* force_divr
    end
end

struct OriginalhPFInteractions <: NonBondInteraction 
    κ::Float64
    N::Vector{Int64}
end

function apply_nonbonds!(args::system,atoms::Vector{atom},forces::Vector{force},energy::Vector{Float64},hPF::OriginalhPFInteractions)
    mesh=init_mesh(args.box,hPF.N[1],hPF.N[2],hPF.N[3])
     for atomii=1:length(atoms)
        cloudincell!(atoms[atomii].pos,mesh)
    end
    DensityatVertex!(mesh)
    Grad_DensVertex!(mesh)
    ParticleFieldInteraction!(atoms,forces,energy,mesh)
end
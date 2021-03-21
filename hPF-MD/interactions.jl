
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
    potentials::Array{Float64,1}
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
    mesh::Mesh
end


function ParticleFieldInteraction!(atoms::Vector{atom},forces::Vector{force},energy::Vector{Float64},mesh::Mesh,hPF::OriginalhPFInteractions)
    avg_den=length(atoms)/prod(mesh.boxsize)
    pic=zeros(8)
    vertex_=zeros(Int64,8)
    for atomii =1:length(atoms)
        ix::Int64=floor(atoms[atomii].pos[1] / mesh.edge_size[1])
        iy::Int64=floor(atoms[atomii].pos[2] / mesh.edge_size[2])
        iz::Int64=floor(atoms[atomii].pos[3] / mesh.edge_size[3])
        index_cell::Int64=getCellIndex([ix,iy,iz],mesh.N)
        cellpos_low::Vector{Float64}=atoms[atomii].pos .- floor.(atoms[atomii].pos ./ mesh.edge_size) .* mesh.edge_size

        pic[1]= (mesh.edge_size[1]-cellpos_low[1]) * (mesh.edge_size[2]-cellpos_low[2]) * (mesh.edge_size[3]-cellpos_low[3]) / prod(mesh.edge_size)
        pic[2]= cellpos_low[1]                     * (mesh.edge_size[2]-cellpos_low[2]) * (mesh.edge_size[3]-cellpos_low[3]) / prod(mesh.edge_size)
        pic[3]= (mesh.edge_size[1]-cellpos_low[1]) * cellpos_low[2]                     * (mesh.edge_size[3]-cellpos_low[3]) / prod(mesh.edge_size)
        pic[4]= cellpos_low[1]                     * cellpos_low[2]                     * (mesh.edge_size[3]-cellpos_low[3]) / prod(mesh.edge_size)
        pic[5]= (mesh.edge_size[1]-cellpos_low[1]) * (mesh.edge_size[2]-cellpos_low[2]) *                     cellpos_low[3] / prod(mesh.edge_size)
        pic[6]= cellpos_low[1]                     * (mesh.edge_size[2]-cellpos_low[2]) *                     cellpos_low[3] / prod(mesh.edge_size)
        pic[7]= (mesh.edge_size[1]-cellpos_low[1]) * cellpos_low[2]                     *                     cellpos_low[3] / prod(mesh.edge_size)
        pic[8]= cellpos_low[1]                     * cellpos_low[2]                     *                     cellpos_low[3] / prod(mesh.edge_size)

        icell=ix
        jcell=iy
        kcell=iz
        icell_plus=ix+1
        jcell_plus=iy+1
        kcell_plus=iz+1

        
        vertex_[1]=getCellIndex(pbc_mesh([icell,jcell,kcell],mesh.N),mesh.N)                 
        vertex_[2]=getCellIndex(pbc_mesh([icell_plus,jcell,kcell],mesh.N),mesh.N)            
        vertex_[3]=getCellIndex(pbc_mesh([icell,jcell_plus,kcell],mesh.N),mesh.N)            
        vertex_[4]=getCellIndex(pbc_mesh([icell_plus,jcell_plus,kcell],mesh.N),mesh.N)       
        vertex_[5]=getCellIndex(pbc_mesh([icell,jcell,kcell_plus],mesh.N),mesh.N)            
        vertex_[6]=getCellIndex(pbc_mesh([icell_plus,jcell,kcell_plus],mesh.N),mesh.N)       
        vertex_[7]=getCellIndex(pbc_mesh([icell,jcell_plus,kcell_plus],mesh.N),mesh.N)       
        vertex_[8]=getCellIndex(pbc_mesh([icell_plus,jcell_plus,kcell_plus],mesh.N),mesh.N)  

        den_atomii=0
        densgrads_atomii=zeros(3)
        for vertexii=1:8
            den_atomii          +=  mesh.vertexes[vertex_[vertexii]]*pic[vertexii]
            densgrads_atomii    +=   mesh.densgrads[vertex_[vertexii],:] .* pic[vertexii]
        end
        energy[atomii] += 1/hPF.κ*(den_atomii-avg_den)/avg_den
        forces[atomii].force += -1 ./ hPF.κ .*densgrads_atomii/avg_den
    end
   
end

function apply_nonbonds!(args::system,atoms::Vector{atom},forces::Vector{force},energy::Vector{Float64},hPF::OriginalhPFInteractions)
    clear_mesh!(hPF.mesh)
    for atomii=1:length(atoms)
        cloudincell!(atoms[atomii].pos,hPF.mesh)
    end
    DensityatVertex!(hPF.mesh)
    Grad_DensVertex!(hPF.mesh)
    ParticleFieldInteraction!(atoms,forces,energy,hPF.mesh,hPF)
end
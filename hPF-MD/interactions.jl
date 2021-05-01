
export TableBondInteractions,OriginalhPFInteractions,FENEBondInteractions,HarmonicBondInteractions


function clear_forces!(forces::Vector{force})
    for forceii=1:length(forces)
        forces[forceii].force[1]=0.0
        forces[forceii].force[2]=0.0
        forces[forceii].force[3]=0.0
    end
end

function clear_energy!(energy::Vector{Float64})
    for i=1:length(energy)
        energy[i]=0.0
    end
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
    # 0.5×k×(r-r0)^2
    for bondii=1:length(bonds)
        delta_d = atoms[bonds[bondii].atom_1+1].pos .- atoms[bonds[bondii].atom_2+1].pos

        # consider pbc
        delta_d .-= args.box .* RInt.( delta_d ./ args.box)

        r=norm(delta_d)
        
        force_=harmonicbond.κ*(harmonicbond.r/r-1.0)
        force_divr = force_/r
        
        # calculate energy
        energy[bonds[bondii].atom_1+1]+=0.5*0.5*harmonicbond.κ*(r-harmonicbond.r)^2
        energy[bonds[bondii].atom_2+1]+=0.5*0.5*harmonicbond.κ*(r-harmonicbond.r)^2

        # calculate forces
        forces[bonds[bondii].atom_1+1].force .+= delta_d .* force_divr
        forces[bonds[bondii].atom_2+1].force .+= -delta_d .* force_divr
    end
end

struct FENEBondInteractions <: BondInteraction 
    K::Float64 #    (energy/distance^2)
    R0::Float64 #   (distance)
    ϵ::Float64 #    (energy)
    σ::Float64 #    (distance)
end

function apply_bonds!(args::system,atoms::Vector{atom},forces::Vector{force},energy::Vector{Float64},bonds::Vector{bond},fenebond::FENEBondInteractions) 

    for bondii=1:length(bonds)
        delta_d = atoms[bonds[bondii].atom_1+1].pos .- atoms[bonds[bondii].atom_2+1].pos

        # consider pbc
        delta_d .-= args.box .* RInt.( delta_d ./ args.box)

        r = (delta_d[1]^2+delta_d[2]^2+delta_d[3]^2)^0.5
        rsq=r^2

        r0sq=fenebond.R0*fenebond.R0
        rlogarg=1-rsq/r0sq

        fbond = -fenebond.K/rlogarg

        if rsq<fenebond.σ^2
            sr2=fenebond.σ^2/rsq
            sr6=sr2*sr2*sr2
            force_divr += 48*fenebond.ϵ*sr6*(sr6-0.5)/rsq
        end

        E=-0.5*fenebond.K*fenebond.R0^2*log(1-(r/fenebond.R0)^2)+4*fenebond.ϵ*((fenebond.σ/r)^12-(fenebond.σ/r)^6)+fenebond.ϵ

        energy[bonds[bondii].atom_1+1]+=E*0.5
        energy[bonds[bondii].atom_2+1]+=E*0.5

        # calculate forces
        forces[bonds[bondii].atom_1+1].force .+= delta_d .* force_divr
        forces[bonds[bondii].atom_2+1].force .+= -delta_d .* force_divr
    end
end

struct OriginalhPFInteractions <: NonBondInteraction 
    κ::Float64
    mesh::Mesh
end

function apply_nonbonds!(args::system,atoms::Vector{atom},forces::Vector{force},energy::Vector{Float64},::NoNonBondInteractions) end

function ParticleFieldInteraction!(atoms::Vector{atom},forces::Vector{force},energy::Vector{Float64},mesh::Mesh,hPF::OriginalhPFInteractions)
    avg_den=length(atoms)/prod(mesh.boxsize)
    pic=zeros(8)
    vertex_=zeros(Int64,8)
    densgrads_atomii=zeros(3)
    for atomii =1:length(atoms)

        ix::Int64=floor(atoms[atomii].pos[1] / mesh.edge_size[1])
        iy::Int64=floor(atoms[atomii].pos[2] / mesh.edge_size[2])
        iz::Int64=floor(atoms[atomii].pos[3] / mesh.edge_size[3])

        index_cell::Int64=getCellIndex(ix,iy,iz,mesh.N)

        δx=atoms[atomii].pos[1] - ix * mesh.edge_size[1]
        δy=atoms[atomii].pos[2] - iy * mesh.edge_size[2]
        δz=atoms[atomii].pos[3] - iz * mesh.edge_size[3]

        pic[1]= (mesh.edge_size[1]-δx) * (mesh.edge_size[2]-δy) * (mesh.edge_size[3]-δz) / prod(mesh.edge_size)
        pic[2]= δx                     * (mesh.edge_size[2]-δy) * (mesh.edge_size[3]-δz) / prod(mesh.edge_size)
        pic[3]= (mesh.edge_size[1]-δx) * δy                     * (mesh.edge_size[3]-δz) / prod(mesh.edge_size)
        pic[4]= δx                     * δy                     * (mesh.edge_size[3]-δz) / prod(mesh.edge_size)
        pic[5]= (mesh.edge_size[1]-δx) * (mesh.edge_size[2]-δy) *                     δz / prod(mesh.edge_size)
        pic[6]= δx                     * (mesh.edge_size[2]-δy) *                     δz / prod(mesh.edge_size)
        pic[7]= (mesh.edge_size[1]-δx) * δy                     *                     δz / prod(mesh.edge_size)
        pic[8]= δx                     * δy                     *                     δz / prod(mesh.edge_size)

        icell=ix
        jcell=iy
        kcell=iz

        icell=pbc_mesh(icell,mesh.N[1])
        jcell=pbc_mesh(jcell,mesh.N[2])
        kcell=pbc_mesh(kcell,mesh.N[3])

        icell_plus=pbc_mesh(icell+1,mesh.N[1])
        jcell_plus=pbc_mesh(jcell+1,mesh.N[1])
        kcell_plus=pbc_mesh(kcell+1,mesh.N[1])

        vertex_[1]=getCellIndex(icell,jcell,kcell,mesh.N)                 
        vertex_[2]=getCellIndex(icell_plus,jcell,kcell,mesh.N)            
        vertex_[3]=getCellIndex(icell,jcell_plus,kcell,mesh.N)            
        vertex_[4]=getCellIndex(icell_plus,jcell_plus,kcell,mesh.N)       
        vertex_[5]=getCellIndex(icell,jcell,kcell_plus,mesh.N)            
        vertex_[6]=getCellIndex(icell_plus,jcell,kcell_plus,mesh.N)       
        vertex_[7]=getCellIndex(icell,jcell_plus,kcell_plus,mesh.N)       
        vertex_[8]=getCellIndex(icell_plus,jcell_plus,kcell_plus,mesh.N)  

        den_atomii=0
        
        densgrads_atomii[1]=0.0
        densgrads_atomii[2]=0.0
        densgrads_atomii[3]=0.0

        for vertexii=1:8
            den_atomii          +=  mesh.vertexes[vertex_[vertexii]]*pic[vertexii]
            densgrads_atomii[1]    +=   mesh.densgrads[vertex_[vertexii],1] * pic[vertexii]
            densgrads_atomii[2]    +=   mesh.densgrads[vertex_[vertexii],2] * pic[vertexii]
            densgrads_atomii[3]    +=   mesh.densgrads[vertex_[vertexii],3] * pic[vertexii]
        end

        energy[atomii] += 1/hPF.κ*(den_atomii-avg_den)/avg_den

        forces[atomii].force[1] += -1 / hPF.κ * densgrads_atomii[1]/avg_den
        forces[atomii].force[2] += -1 / hPF.κ * densgrads_atomii[2]/avg_den
        forces[atomii].force[3] += -1 / hPF.κ * densgrads_atomii[3]/avg_den
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


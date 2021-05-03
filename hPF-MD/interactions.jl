
export TableBondInteractions,FENEBondInteractions,HarmonicBondInteractions,spectralhPFInteractions,originalhPFInteractions


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

struct spectralhPFInteractions <: NonBondInteraction 
    mesh::Mesh
    κ::Float64
    σ::Float64
end

struct originalhPFInteractions <: NonBondInteraction 
    mesh::Mesh
    κ::Float64
end

function apply_nonbonds!(args::system,atoms::Vector{atom},forces::Vector{force},energy::Vector{Float64},::NoNonBondInteractions) end

function ParticleFieldInteraction!(atoms::Vector{atom},forces::Vector{force},energy::Vector{Float64},mesh::Mesh,hPF::originalhPFInteractions)
    avg_den=length(atoms)/prod(mesh.boxsize)
    pic=zeros(8)
    
    # force on the grids: -1 * derivative of the potential on grids
    grid_grad_x,grid_grad_y,grid_grad_z=Grad_DensVertex(mesh)

    for atomii =1:length(atoms)

        ix::Int64=floor(atoms[atomii].pos[1] / mesh.edge_size[1])
        iy::Int64=floor(atoms[atomii].pos[2] / mesh.edge_size[2])
        iz::Int64=floor(atoms[atomii].pos[3] / mesh.edge_size[3])

        δx=atoms[atomii].pos[1] - ix * mesh.edge_size[1]
        δy=atoms[atomii].pos[2] - iy * mesh.edge_size[2]
        δz=atoms[atomii].pos[3] - iz * mesh.edge_size[3]

        pic[1]= (mesh.edge_size[1]-δx) * (mesh.edge_size[2]-δy) * (mesh.edge_size[3]-δz) / prod(mesh.edge_size)
        pic[2]= (mesh.edge_size[1]-δx) * (δy) * (mesh.edge_size[3]-δz) / prod(mesh.edge_size)                 
        pic[3]= (mesh.edge_size[1]-δx) * (δy) * (δz) / prod(mesh.edge_size)                                   
        pic[4]= (mesh.edge_size[1]-δx) * (mesh.edge_size[2]-δy) * (δz) / prod(mesh.edge_size)                 
        pic[5]= (δx) * (mesh.edge_size[2]-δy) * (mesh.edge_size[3]-δz) / prod(mesh.edge_size)                  
        pic[6]= (δx) * (mesh.edge_size[2]-δy) * (δz) / prod(mesh.edge_size)                                   
        pic[7]= (δx) * (δy) * (mesh.edge_size[3]-δz) / prod(mesh.edge_size)                                   
        pic[8]= (δx) * (δy) * (δz) / prod(mesh.edge_size)                                                      

        icell=ix
        jcell=iy
        kcell=iz

        cell_index_x=pbc_mesh(icell,mesh.N[1])
        cell_index_y=pbc_mesh(jcell,mesh.N[2])
        cell_index_z=pbc_mesh(kcell,mesh.N[3])

        cell_index_x_plus=pbc_mesh(cell_index_x+1,mesh.N[1])
        cell_index_y_plus=pbc_mesh(cell_index_y+1,mesh.N[2])
        cell_index_z_plus=pbc_mesh(cell_index_z+1,mesh.N[3])

        cell_index_x=cell_index_x+1
        cell_index_y=cell_index_y+1
        cell_index_z=cell_index_z+1

        cell_index_x_plus=cell_index_x_plus+1
        cell_index_y_plus=cell_index_y_plus+1
        cell_index_z_plus=cell_index_z_plus+1

        energy[atomii]+=mesh.grids[cell_index_x,cell_index_y,cell_index_z]                *pic[1]*1/(avg_den*hPF.κ)
        energy[atomii]+=mesh.grids[cell_index_x,cell_index_y_plus,cell_index_z]           *pic[2]*1/(avg_den*hPF.κ)
        energy[atomii]+=mesh.grids[cell_index_x,cell_index_y_plus,cell_index_z_plus]      *pic[3]*1/(avg_den*hPF.κ)
        energy[atomii]+=mesh.grids[cell_index_x,cell_index_y,cell_index_z_plus]           *pic[4]*1/(avg_den*hPF.κ)
        energy[atomii]+=mesh.grids[cell_index_x_plus,cell_index_y,cell_index_z]           *pic[5]*1/(avg_den*hPF.κ)
        energy[atomii]+=mesh.grids[cell_index_x_plus,cell_index_y,cell_index_z_plus]      *pic[6]*1/(avg_den*hPF.κ)
        energy[atomii]+=mesh.grids[cell_index_x_plus,cell_index_y_plus,cell_index_z]      *pic[7]*1/(avg_den*hPF.κ)
        energy[atomii]+=mesh.grids[cell_index_x_plus,cell_index_y_plus,cell_index_z_plus] *pic[8]*1/(avg_den*hPF.κ)

        energy[atomii]-=1/hPF.κ

        forces[atomii].force[1]+=  grid_grad_x[cell_index_x,cell_index_y,cell_index_z]                *pic[1] * -1 / hPF.κ/avg_den
        forces[atomii].force[1]+=  grid_grad_x[cell_index_x,cell_index_y_plus,cell_index_z]           *pic[2] * -1 / hPF.κ/avg_den
        forces[atomii].force[1]+=  grid_grad_x[cell_index_x,cell_index_y_plus,cell_index_z_plus]      *pic[3] * -1 / hPF.κ/avg_den
        forces[atomii].force[1]+=  grid_grad_x[cell_index_x,cell_index_y,cell_index_z_plus]           *pic[4] * -1 / hPF.κ/avg_den
        forces[atomii].force[1]+=  grid_grad_x[cell_index_x_plus,cell_index_y,cell_index_z]           *pic[5] * -1 / hPF.κ/avg_den
        forces[atomii].force[1]+=  grid_grad_x[cell_index_x_plus,cell_index_y,cell_index_z_plus]      *pic[6] * -1 / hPF.κ/avg_den
        forces[atomii].force[1]+=  grid_grad_x[cell_index_x_plus,cell_index_y_plus,cell_index_z]      *pic[7] * -1 / hPF.κ/avg_den
        forces[atomii].force[1]+=  grid_grad_x[cell_index_x_plus,cell_index_y_plus,cell_index_z_plus] *pic[8] * -1 / hPF.κ/avg_den

        forces[atomii].force[2]+=  grid_grad_y[cell_index_x,cell_index_y,cell_index_z]                *pic[1] * -1 / hPF.κ/avg_den
        forces[atomii].force[2]+=  grid_grad_y[cell_index_x,cell_index_y_plus,cell_index_z]           *pic[2] * -1 / hPF.κ/avg_den
        forces[atomii].force[2]+=  grid_grad_y[cell_index_x,cell_index_y_plus,cell_index_z_plus]      *pic[3] * -1 / hPF.κ/avg_den
        forces[atomii].force[2]+=  grid_grad_y[cell_index_x,cell_index_y,cell_index_z_plus]           *pic[4] * -1 / hPF.κ/avg_den
        forces[atomii].force[2]+=  grid_grad_y[cell_index_x_plus,cell_index_y,cell_index_z]           *pic[5] * -1 / hPF.κ/avg_den
        forces[atomii].force[2]+=  grid_grad_y[cell_index_x_plus,cell_index_y,cell_index_z_plus]      *pic[6] * -1 / hPF.κ/avg_den
        forces[atomii].force[2]+=  grid_grad_y[cell_index_x_plus,cell_index_y_plus,cell_index_z]      *pic[7] * -1 / hPF.κ/avg_den
        forces[atomii].force[2]+=  grid_grad_y[cell_index_x_plus,cell_index_y_plus,cell_index_z_plus] *pic[8] * -1 / hPF.κ/avg_den

        forces[atomii].force[3]+=  grid_grad_z[cell_index_x,cell_index_y,cell_index_z]                *pic[1] * -1 / hPF.κ/avg_den 
        forces[atomii].force[3]+=  grid_grad_z[cell_index_x,cell_index_y_plus,cell_index_z]           *pic[2] * -1 / hPF.κ/avg_den
        forces[atomii].force[3]+=  grid_grad_z[cell_index_x,cell_index_y_plus,cell_index_z_plus]      *pic[3] * -1 / hPF.κ/avg_den
        forces[atomii].force[3]+=  grid_grad_z[cell_index_x,cell_index_y,cell_index_z_plus]           *pic[4] * -1 / hPF.κ/avg_den
        forces[atomii].force[3]+=  grid_grad_z[cell_index_x_plus,cell_index_y,cell_index_z]           *pic[5] * -1 / hPF.κ/avg_den
        forces[atomii].force[3]+=  grid_grad_z[cell_index_x_plus,cell_index_y,cell_index_z_plus]      *pic[6] * -1 / hPF.κ/avg_den
        forces[atomii].force[3]+=  grid_grad_z[cell_index_x_plus,cell_index_y_plus,cell_index_z]      *pic[7] * -1 / hPF.κ/avg_den
        forces[atomii].force[3]+=  grid_grad_z[cell_index_x_plus,cell_index_y_plus,cell_index_z_plus] *pic[8] * -1 / hPF.κ/avg_den

    end
   
end

function ParticleFieldInteraction!(atoms::Vector{atom},forces::Vector{force},energy::Vector{Float64},mesh::Mesh,hPF::spectralhPFInteractions)
    avg_den=length(atoms)/prod(mesh.boxsize)
    pic=zeros(8)

    # apply gaussian filter

    # forward transform
    complex_field=fft(mesh.grids) 
    #convolution with gaussian
    conv_filed=complex_field .* exp.(-0.5 .* hPF.σ^2 .* mesh.Kmag)
    # backward transform
    grids=real(ifft(conv_filed)) 
    
    # force on the grids: -1 * derivative of the potential on grids
    grid_grad_x=real(ifft( conv_filed .* 1/(avg_den*hPF.κ) .* -1 .* mesh.Ls .* 1im)) 
    grid_grad_y=real(ifft( conv_filed .* 1/(avg_den*hPF.κ) .* -1 .* mesh.Ks .* 1im)) 
    grid_grad_z=real(ifft( conv_filed .* 1/(avg_den*hPF.κ) .* -1 .* mesh.Zs .* 1im)) 

    for atomii =1:length(atoms)

        ix::Int64=floor(atoms[atomii].pos[1] / mesh.edge_size[1])
        iy::Int64=floor(atoms[atomii].pos[2] / mesh.edge_size[2])
        iz::Int64=floor(atoms[atomii].pos[3] / mesh.edge_size[3])

        δx=atoms[atomii].pos[1] - ix * mesh.edge_size[1]
        δy=atoms[atomii].pos[2] - iy * mesh.edge_size[2]
        δz=atoms[atomii].pos[3] - iz * mesh.edge_size[3]

        pic[1]= (mesh.edge_size[1]-δx) * (mesh.edge_size[2]-δy) * (mesh.edge_size[3]-δz) / prod(mesh.edge_size)
        pic[2]= (mesh.edge_size[1]-δx) * (δy) * (mesh.edge_size[3]-δz) / prod(mesh.edge_size)                 
        pic[3]= (mesh.edge_size[1]-δx) * (δy) * (δz) / prod(mesh.edge_size)                                   
        pic[4]= (mesh.edge_size[1]-δx) * (mesh.edge_size[2]-δy) * (δz) / prod(mesh.edge_size)                 
        pic[5]= (δx) * (mesh.edge_size[2]-δy) * (mesh.edge_size[3]-δz) / prod(mesh.edge_size)                  
        pic[6]= (δx) * (mesh.edge_size[2]-δy) * (δz) / prod(mesh.edge_size)                                   
        pic[7]= (δx) * (δy) * (mesh.edge_size[3]-δz) / prod(mesh.edge_size)                                   
        pic[8]= (δx) * (δy) * (δz) / prod(mesh.edge_size)                                                      

        icell=ix
        jcell=iy
        kcell=iz

        cell_index_x=pbc_mesh(icell,mesh.N[1])
        cell_index_y=pbc_mesh(jcell,mesh.N[2])
        cell_index_z=pbc_mesh(kcell,mesh.N[3])

        cell_index_x_plus=pbc_mesh(cell_index_x+1,mesh.N[1])
        cell_index_y_plus=pbc_mesh(cell_index_y+1,mesh.N[2])
        cell_index_z_plus=pbc_mesh(cell_index_z+1,mesh.N[3])

        cell_index_x=cell_index_x+1
        cell_index_y=cell_index_y+1
        cell_index_z=cell_index_z+1

        cell_index_x_plus=cell_index_x_plus+1
        cell_index_y_plus=cell_index_y_plus+1
        cell_index_z_plus=cell_index_z_plus+1

        energy[atomii]+=grids[cell_index_x,cell_index_y,cell_index_z]                *pic[1]*1/(avg_den*hPF.κ)
        energy[atomii]+=grids[cell_index_x,cell_index_y_plus,cell_index_z]           *pic[2]*1/(avg_den*hPF.κ)
        energy[atomii]+=grids[cell_index_x,cell_index_y_plus,cell_index_z_plus]      *pic[3]*1/(avg_den*hPF.κ)
        energy[atomii]+=grids[cell_index_x,cell_index_y,cell_index_z_plus]           *pic[4]*1/(avg_den*hPF.κ)
        energy[atomii]+=grids[cell_index_x_plus,cell_index_y,cell_index_z]           *pic[5]*1/(avg_den*hPF.κ)
        energy[atomii]+=grids[cell_index_x_plus,cell_index_y,cell_index_z_plus]      *pic[6]*1/(avg_den*hPF.κ)
        energy[atomii]+=grids[cell_index_x_plus,cell_index_y_plus,cell_index_z]      *pic[7]*1/(avg_den*hPF.κ)
        energy[atomii]+=grids[cell_index_x_plus,cell_index_y_plus,cell_index_z_plus] *pic[8]*1/(avg_den*hPF.κ)
        energy[atomii]-=1/hPF.κ

        forces[atomii].force[1]+=  grid_grad_x[cell_index_x,cell_index_y,cell_index_z]                *pic[1]
        forces[atomii].force[1]+=  grid_grad_x[cell_index_x,cell_index_y_plus,cell_index_z]           *pic[2]
        forces[atomii].force[1]+=  grid_grad_x[cell_index_x,cell_index_y_plus,cell_index_z_plus]      *pic[3]
        forces[atomii].force[1]+=  grid_grad_x[cell_index_x,cell_index_y,cell_index_z_plus]           *pic[4]
        forces[atomii].force[1]+=  grid_grad_x[cell_index_x_plus,cell_index_y,cell_index_z]           *pic[5]
        forces[atomii].force[1]+=  grid_grad_x[cell_index_x_plus,cell_index_y,cell_index_z_plus]      *pic[6]
        forces[atomii].force[1]+=  grid_grad_x[cell_index_x_plus,cell_index_y_plus,cell_index_z]      *pic[7]
        forces[atomii].force[1]+=  grid_grad_x[cell_index_x_plus,cell_index_y_plus,cell_index_z_plus] *pic[8]

        forces[atomii].force[2]+=  grid_grad_y[cell_index_x,cell_index_y,cell_index_z]                *pic[1]
        forces[atomii].force[2]+=  grid_grad_y[cell_index_x,cell_index_y_plus,cell_index_z]           *pic[2]
        forces[atomii].force[2]+=  grid_grad_y[cell_index_x,cell_index_y_plus,cell_index_z_plus]      *pic[3]
        forces[atomii].force[2]+=  grid_grad_y[cell_index_x,cell_index_y,cell_index_z_plus]           *pic[4]
        forces[atomii].force[2]+=  grid_grad_y[cell_index_x_plus,cell_index_y,cell_index_z]           *pic[5]
        forces[atomii].force[2]+=  grid_grad_y[cell_index_x_plus,cell_index_y,cell_index_z_plus]      *pic[6]
        forces[atomii].force[2]+=  grid_grad_y[cell_index_x_plus,cell_index_y_plus,cell_index_z]      *pic[7]
        forces[atomii].force[2]+=  grid_grad_y[cell_index_x_plus,cell_index_y_plus,cell_index_z_plus] *pic[8]

        forces[atomii].force[3]+=  grid_grad_z[cell_index_x,cell_index_y,cell_index_z]                *pic[1]
        forces[atomii].force[3]+=  grid_grad_z[cell_index_x,cell_index_y_plus,cell_index_z]           *pic[2]
        forces[atomii].force[3]+=  grid_grad_z[cell_index_x,cell_index_y_plus,cell_index_z_plus]      *pic[3]
        forces[atomii].force[3]+=  grid_grad_z[cell_index_x,cell_index_y,cell_index_z_plus]           *pic[4]
        forces[atomii].force[3]+=  grid_grad_z[cell_index_x_plus,cell_index_y,cell_index_z]           *pic[5]
        forces[atomii].force[3]+=  grid_grad_z[cell_index_x_plus,cell_index_y,cell_index_z_plus]      *pic[6]
        forces[atomii].force[3]+=  grid_grad_z[cell_index_x_plus,cell_index_y_plus,cell_index_z]      *pic[7]
        forces[atomii].force[3]+=  grid_grad_z[cell_index_x_plus,cell_index_y_plus,cell_index_z_plus] *pic[8]

    end
   
end


function apply_nonbonds!(args::system,atoms::Vector{atom},forces::Vector{force},energy::Vector{Float64},hPF::NonBondInteraction)

    clear_mesh!(hPF.mesh)
    for atomii=1:length(atoms)
        cloudincell!(atoms[atomii].pos,hPF.mesh)
    end
    ParticleFieldInteraction!(atoms,forces,energy,hPF.mesh,hPF)
end

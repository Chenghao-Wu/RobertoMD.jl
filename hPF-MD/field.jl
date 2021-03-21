
export Mesh,clear_mesh!,init_mesh

struct Mesh
    edge_size::Array{Float64,1}
    N::Array{Int64,1}
    cells::Array{Float64,2}
    vertexes::Array{Float64,1}
    densgrads::Array{Float64,2}
    boxsize::Array{Float64,1}
end

function init_mesh(boxsize::Vector{Float64},Nx::Int64,Ny::Int64,Nz::Int64)
    egde_x=boxsize[1]/Nx
    egde_y=boxsize[2]/Ny
    egde_z=boxsize[3]/Nz
    edge_size=[egde_x,egde_y,egde_z]
    N=[Nx,Ny,Nz]
    N_total::Int64=Nx*Ny*Nz
    cells=zeros(N_total,8)
    vertexes=zeros(N_total) 
    densgrads=zeros(N_total,3) 
    return Mesh(edge_size,N,cells,vertexes,densgrads,boxsize)
end

function clear_mesh!(mesh::Mesh)
    for i=1:prod(mesh.N)
        mesh.cells[i,:]=zeros(1,8)
        mesh.vertexes[i]=0.0
        mesh.densgrads[i,:]=zeros(1,3) 
    end
end

function getCellIndex(ixyz::Vector{Int64},N::Vector{Int64})
    cellindex::Int64=1+ixyz[1]+ixyz[2]*N[1]+ixyz[3]*N[1]*N[2]
    return cellindex
end

function getiXYZfromCellIndex(cellindex::Int64,N::Vector{Int64})
    index_z::Int64=floor((cellindex-1)/(N[1]*N[2]))
    index_y::Int64=floor(((cellindex-1)-N[1]*N[2]*index_z)/N[1])
    index_x::Int64=((cellindex-1)-N[1]*N[2]*index_z)%N[1]
    return [index_x,index_y,index_z]
end

function pbc_mesh(index::Vector{Int64},N_mesh::Vector{Int64})
    if index[1] > N_mesh[1]-1
        index[1]-=N_mesh[1]-1
    end
    if index[2] > N_mesh[2]-1
        index[2]-=N_mesh[2]-1
    end
    if index[3] > N_mesh[3]-1
        index[3]-=N_mesh[3]-1
    end
    if index[1] < 0
        index[1]+=N_mesh[1]-1
    end
    if index[2] < 0
        index[2]+=N_mesh[2]-1
    end
    if index[3] < 0
        index[3]+=N_mesh[3]-1
    end
    return index
end

function pbc_particle(position::Vector{Float64},mesh::Mesh)
    position .-= mesh.boxsize .* floor.( position ./ mesh.boxsize)
    return position
end

function cloudincell!(position::Vector{Float64},mesh::Mesh)

    # cell index = 1 + floor(x/cell_sizex) + floor(y/cell_sizey)*Nx + floor(z/cell_sizez)*Nx*Ny 

    # 3           4
    # -------------
    # |  |  |  |  | 
    # |  |  |  |  | 
    # |  |  |  |  |
    # -------------
    # 1           2 
    
    position=pbc_particle(position,mesh)

    ix,iy,iz=floor.(Int64, position ./ mesh.edge_size)

    index_cell::Int64=getCellIndex([ix,iy,iz],mesh.N)
    cellpos_low::Vector{Float64}=position .- floor.(position ./ mesh.edge_size) .* mesh.edge_size
    #delta_pos=position .- index_cell.*mesh.edge_size
    gridsize=mesh.edge_size[1]*mesh.edge_size[2]*mesh.edge_size[3]
    
    mesh.cells[index_cell,1] += (mesh.edge_size[1]-cellpos_low[1]) * (mesh.edge_size[2]-cellpos_low[2]) * (mesh.edge_size[3]-cellpos_low[3]) / gridsize
    mesh.cells[index_cell,2] += cellpos_low[1]                     * (mesh.edge_size[2]-cellpos_low[2]) * (mesh.edge_size[3]-cellpos_low[3]) / gridsize
    mesh.cells[index_cell,3] += (mesh.edge_size[1]-cellpos_low[1]) * cellpos_low[2]                     * (mesh.edge_size[3]-cellpos_low[3]) / gridsize
    mesh.cells[index_cell,4] += cellpos_low[1]                     * cellpos_low[2]                     * (mesh.edge_size[3]-cellpos_low[3]) / gridsize
    mesh.cells[index_cell,5] += (mesh.edge_size[1]-cellpos_low[1]) * (mesh.edge_size[2]-cellpos_low[2]) *                     cellpos_low[3] / gridsize
    mesh.cells[index_cell,6] += cellpos_low[1]                     * (mesh.edge_size[2]-cellpos_low[2]) *                     cellpos_low[3] / gridsize
    mesh.cells[index_cell,7] += (mesh.edge_size[1]-cellpos_low[1]) * cellpos_low[2]                     *                     cellpos_low[3] / gridsize
    mesh.cells[index_cell,8] += cellpos_low[1]                     * cellpos_low[2]                     *                     cellpos_low[3] / gridsize
    
end

function move_front(N::Vector{Int64})
    return 0
end

function DensityatVertex!(mesh::Mesh)
    for cellii=1:prod(mesh.N)
        icell,jcell,kcell=getiXYZfromCellIndex(cellii,mesh.N)
        icell_plus=icell+1
        jcell_plus=jcell+1
        kcell_plus=kcell+1

        mesh.vertexes[getCellIndex(pbc_mesh([icell,jcell,kcell],mesh.N),mesh.N)]                 +=      mesh.cells[cellii,1]
        mesh.vertexes[getCellIndex(pbc_mesh([icell_plus,jcell,kcell],mesh.N),mesh.N)]            +=      mesh.cells[cellii,2]
        mesh.vertexes[getCellIndex(pbc_mesh([icell,jcell_plus,kcell],mesh.N),mesh.N)]            +=      mesh.cells[cellii,3]
        mesh.vertexes[getCellIndex(pbc_mesh([icell_plus,jcell_plus,kcell],mesh.N),mesh.N)]       +=      mesh.cells[cellii,4]
        mesh.vertexes[getCellIndex(pbc_mesh([icell,jcell,kcell_plus],mesh.N),mesh.N)]            +=      mesh.cells[cellii,5]
        mesh.vertexes[getCellIndex(pbc_mesh([icell_plus,jcell,kcell_plus],mesh.N),mesh.N)]       +=      mesh.cells[cellii,6]
        mesh.vertexes[getCellIndex(pbc_mesh([icell,jcell_plus,kcell_plus],mesh.N),mesh.N)]       +=      mesh.cells[cellii,7]
        mesh.vertexes[getCellIndex(pbc_mesh([icell_plus,jcell_plus,kcell_plus],mesh.N),mesh.N)]  +=      mesh.cells[cellii,8]
    end
    
end

function Grad_DensVertex!(mesh::Mesh)
    for cellii=1:prod(mesh.N)
        icell,jcell,kcell=getiXYZfromCellIndex(cellii,mesh.N)

        icell_plus=icell+1
        jcell_plus=jcell+1
        kcell_plus=kcell+1

        icell_minus=icell-1
        jcell_minus=jcell-1
        kcell_minus=kcell-1

        dens_ii         = getCellIndex(pbc_mesh([icell,jcell,kcell],mesh.N),mesh.N)
        dens_ii_right   = getCellIndex(pbc_mesh([icell_plus,jcell,kcell],mesh.N),mesh.N)
        dens_ii_left    = getCellIndex(pbc_mesh([icell_minus,jcell,kcell],mesh.N),mesh.N)
        dens_ii_front   = getCellIndex(pbc_mesh([icell,jcell_plus,kcell],mesh.N),mesh.N)
        dens_ii_back    = getCellIndex(pbc_mesh([icell,jcell_minus,kcell],mesh.N),mesh.N)
        dens_ii_up      = getCellIndex(pbc_mesh([icell,jcell,kcell_plus],mesh.N),mesh.N)
        dens_ii_down    = getCellIndex(pbc_mesh([icell,jcell,kcell_minus],mesh.N),mesh.N)

        
        grad_x=0.5*(mesh.vertexes[dens_ii_right]-mesh.vertexes[dens_ii_left])
        grad_y=0.5*(mesh.vertexes[dens_ii_front]-mesh.vertexes[dens_ii_back])
        grad_z=0.5*(mesh.vertexes[dens_ii_up]-mesh.vertexes[dens_ii_down])
        mesh.densgrads[dens_ii,:] += [grad_x,grad_y,grad_z]
    end
    
end

function gaussian_filter(σ::Float64,pos::Vector{Float64},pos2::Float64)
    Δr ::Vector{Float64} = pos .- pos2
    gauss::Vector{Float64}=1 ./ (sqrt(2 * π) * σ) .* exp.(-1 .* (Δr ./ 2*σ).^2)
    return gauss
end
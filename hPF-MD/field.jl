
export Mesh,clear_mesh!,init_mesh, meshgrid

struct Mesh
    edge_size::Array{Float64,1}
    N::Array{Int64,1}
    grids::Array{Float64,3}
    boxsize::Array{Float64,1}
    boxmatrix::Array{Float64,2}
    Kmag::Array{Float64,3}
    Ks::Array{Float64,3}
    Ls::Array{Float64,3}
    Zs::Array{Float64,3}
end

function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T},
    vz::AbstractVector{T}) where T
    m, n, o = length(vy), length(vx), length(vz)
    vx = reshape(vx, 1, n, 1)
    vy = reshape(vy, m, 1, 1)
    vz = reshape(vz, 1, 1, o)
    om = ones(Int, m)
    on = ones(Int, n)
    oo = ones(Int, o)
    (vx[om, :, oo], vy[:, on, oo], vz[om, on, :])
end

function init_mesh(boxsize::Vector{Float64},Nx::Int64,Ny::Int64,Nz::Int64)
    egde_x=boxsize[1]/Nx
    egde_y=boxsize[2]/Ny
    egde_z=boxsize[3]/Nz
    edge_size=[egde_x,egde_y,egde_z]
    N=[Nx,Ny,Nz]
    N_total::Int64=Nx*Ny*Nz
    grids=zeros(Nx,Ny,Nz)
    boxmatrix=[boxsize[1] 0 0;0 boxsize[2] 0;0 0 boxsize[3]]

    ks=FFTW.fftfreq(N[1], 2π/(boxsize[1]/N[1]))
    ls=FFTW.fftfreq(N[2], 2π/(boxsize[2]/N[1]))
    zs=FFTW.fftfreq(N[3], 2π/(boxsize[3]/N[1]))

    Ks,Ls,Zs=meshgrid(ks,ls,zs)
    #print(Ks)
    
    Kmag=(Ks.^2 .+Ls.^2 .+ Zs.^2)

    return Mesh(edge_size,N,grids,boxsize,boxmatrix,Kmag,Ks,Ls,Zs)
end

function clear_mesh!(mesh::Mesh)
    for i=1:mesh.N[1]
        for j=1:mesh.N[2]
            for k=1:mesh.N[3]
                mesh.grids[i,j,k]=0.0
                mesh.grids[i,j,k]=0.0
                mesh.grids[i,j,k]=0.0
                mesh.grids[i,j,k]=0.0
                mesh.grids[i,j,k]=0.0
                mesh.grids[i,j,k]=0.0
                mesh.grids[i,j,k]=0.0
                mesh.grids[i,j,k]=0.0
            end
        end
    end
end

function getCellIndex(ix::Int64,iy::Int64,iz::Int64,N::Vector{Int64})
    cellindex::Int64=1+ix+iy*N[1]+iz*N[1]*N[2]
    return cellindex
end

function getiXYZfromCellIndex(ixyz::Vector{Int64},cellindex::Int64,N::Vector{Int64})
    ixyz[3]=floor((cellindex-1)/(N[1]*N[2]))
    ixyz[2]=floor(((cellindex-1)-N[1]*N[2]*ixyz[3])/N[1])
    ixyz[1]=((cellindex-1)-N[1]*N[2]*ixyz[3])%N[1]
    return ixyz
end

function pbc_mesh(x_index::Int64,N_mesh::Int64)
    if x_index > N_mesh-1
        x_index-=N_mesh
        if x_index>N_mesh-1
            x_index-=N_mesh
        end
    end
    if x_index < 0
        x_index+=N_mesh
    end
    return x_index
end

function pbc_particle(position::Vector{Float64},mesh::Mesh)
    position .-= mesh.boxsize .* floor.( position ./ mesh.boxsize)
    return position
end

function pbc_particle(position::Float64,boxsize::Float64)
    position -= boxsize * floor( position / boxsize)
    return position
end

function cloudincell(wrapped_pos::Array{Float64,1},gridedge::Array{Float64,1},meshnumber::Array{Int64,1})
    
    δx_lower=wrapped_pos[1] - floor(wrapped_pos[1] / edge_size[1]) * edge_size[1]
    δy_lower=wrapped_pos[2] - floor(wrapped_pos[2] / edge_size[2]) * edge_size[2]
    δz_lower=wrapped_pos[3] - floor(wrapped_pos[3] / edge_size[3]) * edge_size[3]

    δx_upper=edge_size[1]-δx_lower
    δy_upper=edge_size[1]-δy_lower
    δz_upper=edge_size[1]-δz_lower

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
    
    position_x=pbc_particle(position[1],mesh.boxsize[1])
    position_y=pbc_particle(position[2],mesh.boxsize[2])
    position_z=pbc_particle(position[3],mesh.boxsize[3])

    δx=position[1] - floor(position[1] / mesh.edge_size[1]) * mesh.edge_size[1]
    δy=position[2] - floor(position[2] / mesh.edge_size[2]) * mesh.edge_size[2]
    δz=position[3] - floor(position[3] / mesh.edge_size[3]) * mesh.edge_size[3]

    cell_index_x=floor(Int64, position_x / mesh.edge_size[1])
    cell_index_y=floor(Int64, position_y / mesh.edge_size[2])
    cell_index_z=floor(Int64, position_z / mesh.edge_size[3])

    cell_index_x_plus=pbc_mesh(cell_index_x+1,mesh.N[1])
    cell_index_y_plus=pbc_mesh(cell_index_y+1,mesh.N[2])
    cell_index_z_plus=pbc_mesh(cell_index_z+1,mesh.N[3])

    cell_index_x=cell_index_x+1
    cell_index_y=cell_index_y+1
    cell_index_z=cell_index_z+1

    cell_index_x_plus=cell_index_x_plus+1
    cell_index_y_plus=cell_index_y_plus+1
    cell_index_z_plus=cell_index_z_plus+1

    mesh.grids[cell_index_x,cell_index_y,cell_index_z]                += (mesh.edge_size[1]-δx) * (mesh.edge_size[2]-δy) * (mesh.edge_size[3]-δz) / prod(mesh.edge_size)
    mesh.grids[cell_index_x,cell_index_y_plus,cell_index_z]           += (mesh.edge_size[1]-δx) * (δy) * (mesh.edge_size[3]-δz) / prod(mesh.edge_size)                 
    mesh.grids[cell_index_x,cell_index_y_plus,cell_index_z_plus]      += (mesh.edge_size[1]-δx) * (δy) * (δz) / prod(mesh.edge_size)                                   
    mesh.grids[cell_index_x,cell_index_y,cell_index_z_plus]           += (mesh.edge_size[1]-δx) * (mesh.edge_size[2]-δy) * (δz) / prod(mesh.edge_size)                 
    mesh.grids[cell_index_x_plus,cell_index_y,cell_index_z]           += (δx) * (mesh.edge_size[2]-δy) * (mesh.edge_size[3]-δz) / prod(mesh.edge_size)                  
    mesh.grids[cell_index_x_plus,cell_index_y,cell_index_z_plus]      += (δx) * (mesh.edge_size[2]-δy) * (δz) / prod(mesh.edge_size)                                   
    mesh.grids[cell_index_x_plus,cell_index_y_plus,cell_index_z]      += (δx) * (δy) * (mesh.edge_size[3]-δz) / prod(mesh.edge_size)                                   
    mesh.grids[cell_index_x_plus,cell_index_y_plus,cell_index_z_plus] += (δx) * (δy) * (δz) / prod(mesh.edge_size)                                                      
    
end

function Grad_DensVertex(mesh::Mesh)
    grid_grad_x=zeros(mesh.N[1],mesh.N[2],mesh.N[3])
    grid_grad_y=zeros(mesh.N[1],mesh.N[2],mesh.N[3])
    grid_grad_z=zeros(mesh.N[1],mesh.N[2],mesh.N[3])

    for icell=0:mesh.N[1]-1
        for jcell=0:mesh.N[2]-1
            for kcell=0:mesh.N[3]-1


                icell_plus=pbc_mesh(icell+1,mesh.N[1])+1
                jcell_plus=pbc_mesh(jcell+1,mesh.N[2])+1
                kcell_plus=pbc_mesh(kcell+1,mesh.N[3])+1

                icell_minus=pbc_mesh(icell-1,mesh.N[1])+1
                jcell_minus=pbc_mesh(jcell-1,mesh.N[2])+1
                kcell_minus=pbc_mesh(kcell-1,mesh.N[3])+1

                icell_=icell+1
                jcell_=jcell+1
                kcell_=kcell+1
                
                grad_x=0.5*(mesh.grids[icell_plus,jcell_,kcell_]-mesh.grids[icell_minus,jcell_,kcell_])
                grad_y=0.5*(mesh.grids[icell_,jcell_plus,kcell_]-mesh.grids[icell_,jcell_minus,kcell_])
                grad_z=0.5*(mesh.grids[icell_,jcell_,kcell_plus]-mesh.grids[icell_,jcell_,kcell_minus])
                grid_grad_x[icell_,jcell_,kcell_] += grad_x
                grid_grad_y[icell_,jcell_,kcell_] += grad_y
                grid_grad_z[icell_,jcell_,kcell_] += grad_z
            end
        end
    end
    return grid_grad_x,grid_grad_y,grid_grad_z
end
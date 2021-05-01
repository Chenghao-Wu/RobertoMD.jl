
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
        mesh.cells[i,1]=0.0
        mesh.cells[i,2]=0.0
        mesh.cells[i,3]=0.0
        mesh.cells[i,4]=0.0
        mesh.cells[i,5]=0.0
        mesh.cells[i,6]=0.0
        mesh.cells[i,7]=0.0
        mesh.cells[i,8]=0.0
        mesh.vertexes[i]=0.0
        mesh.densgrads[i,1]=0.0
        mesh.densgrads[i,2]=0.0
        mesh.densgrads[i,3]=0.0
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

    ix=floor(Int64, position[1] / mesh.edge_size[1])
    iy=floor(Int64, position[2] / mesh.edge_size[2])
    iz=floor(Int64, position[3] / mesh.edge_size[3])

    index_cell::Int64=getCellIndex(ix,iy,iz,mesh.N)
    
    δx=position[1] - floor(position[1] / mesh.edge_size[1]) * mesh.edge_size[1]
    δy=position[2] - floor(position[2] / mesh.edge_size[2]) * mesh.edge_size[2]
    δz=position[3] - floor(position[3] / mesh.edge_size[3]) * mesh.edge_size[3]

    mesh.cells[index_cell,1] += (mesh.edge_size[1]-δx) * (mesh.edge_size[2]-δy) * (mesh.edge_size[3]-δz) / prod(mesh.edge_size)
    mesh.cells[index_cell,2] += δx                     * (mesh.edge_size[2]-δy) * (mesh.edge_size[3]-δz) / prod(mesh.edge_size)
    mesh.cells[index_cell,3] += (mesh.edge_size[1]-δx) * δy                     * (mesh.edge_size[3]-δz) / prod(mesh.edge_size)
    mesh.cells[index_cell,4] += δx                     * δy                     * (mesh.edge_size[3]-δz) / prod(mesh.edge_size)
    mesh.cells[index_cell,5] += (mesh.edge_size[1]-δx) * (mesh.edge_size[2]-δy) *                     δz / prod(mesh.edge_size)
    mesh.cells[index_cell,6] += δx                     * (mesh.edge_size[2]-δy) *                     δz / prod(mesh.edge_size)
    mesh.cells[index_cell,7] += (mesh.edge_size[1]-δx) * δy                     *                     δz / prod(mesh.edge_size)
    mesh.cells[index_cell,8] += δx                     * δy                     *                     δz / prod(mesh.edge_size)
    
end

function gaussian_kernel_3d(sigma=1.0)
    kernel_3D=zeros(3,3,3)
    center=2
    for i=1:3
        for j=1:3
            for z=1:3
                kernel_3D[i,j,z]=1 / (sqrt(2 * π) * sigma)^3 * exp(-(((i - center)^2+(j - center)^2+(z - center)^2 )/ sigma^2/2))
            end
        end
    end
    kernel_3D.*=1/sum(kernel_3D)
    return kernel_3D
end

function gaussian_filter(mesh::Mesh,sigma::Float64)
    gaussian_kernel=gaussian_kernel_3d(sigma)
    new_densityatvertex=zeros(mesh.N[1]*mesh.N[2]*mesh.N[3]) 
    for icell=0:mesh.N[1]-1
        for jcell=0:mesh.N[2]-1
            for kcell=0:mesh.N[3]-1
                filter=gaussian_kernel*mesh.vertexes[getCellIndex(icell,jcell,kcell,mesh.N)]  

                icell_plus=pbc_mesh(icell+1,mesh.N[1])
                jcell_plus=pbc_mesh(jcell+1,mesh.N[2])
                kcell_plus=pbc_mesh(kcell+1,mesh.N[3])

                icell_minus=pbc_mesh(icell-1,mesh.N[1])
                jcell_minus=pbc_mesh(jcell-1,mesh.N[2])
                kcell_minus=pbc_mesh(kcell-1,mesh.N[3])

                new_densityatvertex[getCellIndex(icell_minus,jcell_minus,kcell_minus,mesh.N)]+=filter[1,1,1]
                new_densityatvertex[getCellIndex(icell_minus,jcell_minus,kcell,mesh.N)]+=filter[1,1,2]
                new_densityatvertex[getCellIndex(icell_minus,jcell_minus,kcell_plus,mesh.N)]+=filter[1,1,3]

                new_densityatvertex[getCellIndex(icell_minus,jcell,kcell_minus,mesh.N)]+=filter[1,2,1]
                new_densityatvertex[getCellIndex(icell_minus,jcell,kcell,mesh.N)]+=filter[1,2,2]
                new_densityatvertex[getCellIndex(icell_minus,jcell,kcell_plus,mesh.N)]+=filter[1,2,3]
                
                new_densityatvertex[getCellIndex(icell_minus,jcell_plus,kcell_minus,mesh.N)]+=filter[1,3,1]
                new_densityatvertex[getCellIndex(icell_minus,jcell_plus,kcell,mesh.N)]+=filter[1,3,2]
                new_densityatvertex[getCellIndex(icell_minus,jcell_plus,kcell_plus,mesh.N)]+=filter[1,3,3]

                new_densityatvertex[getCellIndex(icell,jcell_minus,kcell_minus,mesh.N)]+=filter[2,1,1]
                new_densityatvertex[getCellIndex(icell,jcell_minus,kcell,mesh.N)]+=filter[2,1,2]
                new_densityatvertex[getCellIndex(icell,jcell_minus,kcell_plus,mesh.N)]+=filter[2,1,3]

                new_densityatvertex[getCellIndex(icell,jcell,kcell_minus,mesh.N)]+=filter[2,2,1]
                new_densityatvertex[getCellIndex(icell,jcell,kcell,mesh.N)]+=filter[2,2,2]
                new_densityatvertex[getCellIndex(icell,jcell,kcell_plus,mesh.N)]+=filter[2,2,3]

                new_densityatvertex[getCellIndex(icell,jcell_plus,kcell_minus,mesh.N)]+=filter[2,3,1]
                new_densityatvertex[getCellIndex(icell,jcell_plus,kcell,mesh.N)]+=filter[2,3,2]
                new_densityatvertex[getCellIndex(icell,jcell_plus,kcell_plus,mesh.N)]+=filter[2,3,3]

                new_densityatvertex[getCellIndex(icell_plus,jcell_minus,kcell_minus,mesh.N)]+=filter[3,1,1]
                new_densityatvertex[getCellIndex(icell_plus,jcell_minus,kcell,mesh.N)]+=filter[3,1,2]
                new_densityatvertex[getCellIndex(icell_plus,jcell_minus,kcell_plus,mesh.N)]+=filter[3,1,3]

                new_densityatvertex[getCellIndex(icell_plus,jcell,kcell_minus,mesh.N)]+=filter[3,2,1]
                new_densityatvertex[getCellIndex(icell_plus,jcell,kcell,mesh.N)]+=filter[3,2,2]
                new_densityatvertex[getCellIndex(icell_plus,jcell,kcell_plus,mesh.N)]+=filter[3,2,3]

                new_densityatvertex[getCellIndex(icell_plus,jcell_plus,kcell_minus,mesh.N)]+=filter[3,3,1]
                new_densityatvertex[getCellIndex(icell_plus,jcell_plus,kcell,mesh.N)]+=filter[3,3,2]
                new_densityatvertex[getCellIndex(icell_plus,jcell_plus,kcell_plus,mesh.N)]+=filter[3,3,3]
            end
        end
    end
    return new_densityatvertex
end 

function DensityatVertex!(mesh::Mesh)
    ixyz=zeros(Int64,3)
    for cellii=1:prod(mesh.N)

        icell,jcell,kcell=getiXYZfromCellIndex(ixyz,cellii,mesh.N)

        icell=pbc_mesh(icell,mesh.N[1])
        jcell=pbc_mesh(jcell,mesh.N[2])
        kcell=pbc_mesh(kcell,mesh.N[3])

        icell_plus=pbc_mesh(icell+1,mesh.N[1])
        jcell_plus=pbc_mesh(jcell+1,mesh.N[2])
        kcell_plus=pbc_mesh(kcell+1,mesh.N[3])

        mesh.vertexes[getCellIndex(icell,jcell,kcell,mesh.N)]                 +=      mesh.cells[cellii,1]
        mesh.vertexes[getCellIndex(icell_plus,jcell,kcell,mesh.N)]            +=      mesh.cells[cellii,2]
        mesh.vertexes[getCellIndex(icell,jcell_plus,kcell,mesh.N)]            +=      mesh.cells[cellii,3]
        mesh.vertexes[getCellIndex(icell_plus,jcell_plus,kcell,mesh.N)]       +=      mesh.cells[cellii,4]
        mesh.vertexes[getCellIndex(icell,jcell,kcell_plus,mesh.N)]            +=      mesh.cells[cellii,5]
        mesh.vertexes[getCellIndex(icell_plus,jcell,kcell_plus,mesh.N)]       +=      mesh.cells[cellii,6]
        mesh.vertexes[getCellIndex(icell,jcell_plus,kcell_plus,mesh.N)]       +=      mesh.cells[cellii,7]
        mesh.vertexes[getCellIndex(icell_plus,jcell_plus,kcell_plus,mesh.N)]  +=      mesh.cells[cellii,8]
    end

    # apply a 3-d gaussian filter
    #mesh.vertexes=gaussian_filter(mesh,0.5)

end

function Grad_DensVertex!(mesh::Mesh)
    ixyz=zeros(Int64,3)
    for cellii=1:prod(mesh.N)

        icell,jcell,kcell=getiXYZfromCellIndex(ixyz,cellii,mesh.N)

        icell_plus=pbc_mesh(icell+1,mesh.N[1])
        jcell_plus=pbc_mesh(jcell+1,mesh.N[2])
        kcell_plus=pbc_mesh(kcell+1,mesh.N[3])

        icell_minus=pbc_mesh(icell-1,mesh.N[1])
        jcell_minus=pbc_mesh(jcell-1,mesh.N[2])
        kcell_minus=pbc_mesh(kcell-1,mesh.N[3])

        dens_ii         = getCellIndex(icell,jcell,kcell,mesh.N)
        dens_ii_right   = getCellIndex(icell_plus,jcell,kcell,mesh.N)
        dens_ii_left    = getCellIndex(icell_minus,jcell,kcell,mesh.N)
        dens_ii_front   = getCellIndex(icell,jcell_plus,kcell,mesh.N)
        dens_ii_back    = getCellIndex(icell,jcell_minus,kcell,mesh.N)
        dens_ii_up      = getCellIndex(icell,jcell,kcell_plus,mesh.N)
        dens_ii_down    = getCellIndex(icell,jcell,kcell_minus,mesh.N)

        
        grad_x=0.5*(mesh.vertexes[dens_ii_right]-mesh.vertexes[dens_ii_left])
        grad_y=0.5*(mesh.vertexes[dens_ii_front]-mesh.vertexes[dens_ii_back])
        grad_z=0.5*(mesh.vertexes[dens_ii_up]-mesh.vertexes[dens_ii_down])
        mesh.densgrads[dens_ii,1] += grad_x
        mesh.densgrads[dens_ii,2] += grad_y
        mesh.densgrads[dens_ii,3] += grad_z
    end
    
end



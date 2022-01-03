
function force_field(sys::System,field::NoField,comm::MPI.Comm)
    forces=zeros(sys.local_N,3)
    return forces
end

function clear_field!(sys::System)
    for i in 1:sys.num_atomtype
        fill!(sys.local_field[i] , 0.0)    
        fill!(sys.global_field[i] , 0.0)  
        fill!(sys.field_gradient[i] , 0.0)  
    end
end

function SerializeCellIndex(i::Int64,j::Int64,k::Int64,celllength::Int64)
    return 1+(i-1)+(j-1)*celllength+(k-1)*celllength*celllength
end 

function SerializeCellGradIndex(i::Int64,j::Int64,k::Int64,axis::Int64,celllength::Int64)
    return 1+(i-1)+(j-1)*celllength+(k-1)*celllength*celllength+(axis-1)*celllength*celllength*celllength
end 

function particle2cell_vertex(cell_length::Float64,delta_r::Array{Float64,1},index::Int64)
    if index==1
        return (cell_length-delta_r[1]) * (cell_length-delta_r[2]) * (cell_length-delta_r[3])
    elseif index==2
        return (cell_length-delta_r[1]) * (delta_r[2]) * (cell_length-delta_r[3])
    elseif index==3
        return (cell_length-delta_r[1]) * (delta_r[2]) * (delta_r[3])
    elseif index==4
        return (cell_length-delta_r[1]) * (cell_length-delta_r[2]) * (delta_r[3])
    elseif index==5
        return (delta_r[1]) * (cell_length-delta_r[2]) * (cell_length-delta_r[3])
    elseif index ==6
        return (delta_r[1]) * (cell_length-delta_r[2]) * (delta_r[3])
    elseif index==7
        return (delta_r[1]) * (delta_r[2]) * (cell_length-delta_r[3])
    elseif index ==8
        return (delta_r[1]) * (delta_r[2]) * (delta_r[3])
    end
end

function particle2cell!(sys::System)
    
    V_cell=sys.mesh.cell_length*sys.mesh.cell_length*sys.mesh.cell_length
    for i in 1:sys.local_N
        delta_r=sys.coords[i,1:3] .- floor.(sys.coords[i,1:3]/sys.mesh.cell_length) .* sys.mesh.cell_length
        sys.particle2cell[i,1]=particle2cell_vertex(sys.mesh.cell_length,delta_r,1)     
        sys.particle2cell[i,2]=particle2cell_vertex(sys.mesh.cell_length,delta_r,2)             
        sys.particle2cell[i,3]=particle2cell_vertex(sys.mesh.cell_length,delta_r,3)          
        sys.particle2cell[i,4]=particle2cell_vertex(sys.mesh.cell_length,delta_r,4)             
        sys.particle2cell[i,5]=particle2cell_vertex(sys.mesh.cell_length,delta_r,5)              
        sys.particle2cell[i,6]=particle2cell_vertex(sys.mesh.cell_length,delta_r,6)          
        sys.particle2cell[i,7]=particle2cell_vertex(sys.mesh.cell_length,delta_r,7)          
        sys.particle2cell[i,8]=particle2cell_vertex(sys.mesh.cell_length,delta_r,8)
    end
    sys.particle2cell=sys.particle2cell./V_cell
end

function particle2density!(sys::System)

    for i in 1:sys.local_N
        for j in 1:8
            sys.local_field[sys.types[i]][SerializeCellIndex(sys.cellindex[i][j,1],sys.cellindex[i][j,2],sys.cellindex[i][j,3],sys.mesh.num_cells)] += sys.particle2cell[i,j]
        end
    end 
end

function pbc_mesh(index::Int64,num_cells::Int64)
    if index>0
        return index-floor(Int64,index/(num_cells-1))*(num_cells-1)
    else
        return index-floor(Int64,index/(num_cells-1))*(num_cells-1)
    end
end

function update_cellindex!(sys::System)

    cell_index=floor.(Int64, sys.coords / sys.mesh.cell_length) # 0 to num_cells-1, PBC: 0=num_cells-1
    cell_index=cell_index-floor.(Int64,cell_index./((sys.mesh.num_cells-1)))*(sys.mesh.num_cells-1)
    cell_index_plus=cell_index.+1-floor.(Int64,(cell_index.+1)./((sys.mesh.num_cells-1)))*(sys.mesh.num_cells-1)

    cell_index=cell_index.+1
    cell_index_plus=cell_index_plus.+1

    for i in 1:sys.local_N

        sys.cellindex[i][1,1]=cell_index[i,1]
        sys.cellindex[i][1,2]=cell_index[i,2]
        sys.cellindex[i][1,3]=cell_index[i,3]      

        sys.cellindex[i][2,1]=cell_index[i,1]   
        sys.cellindex[i][2,2]=cell_index_plus[i,2]
        sys.cellindex[i][2,3]=cell_index[i,3]           

        sys.cellindex[i][3,1]=cell_index[i,1]  
        sys.cellindex[i][3,2]=cell_index_plus[i,2] 
        sys.cellindex[i][3,3]=cell_index_plus[i,3]

        sys.cellindex[i][4,1]=cell_index[i,1] 
        sys.cellindex[i][4,2]=cell_index[i,2]
        sys.cellindex[i][4,3]=cell_index_plus[i,3]

        sys.cellindex[i][5,1]=cell_index_plus[i,1] 
        sys.cellindex[i][5,2]=cell_index[i,2]  
        sys.cellindex[i][5,3]=cell_index[i,3]   

        sys.cellindex[i][6,1]=cell_index_plus[i,1]
        sys.cellindex[i][6,2]=cell_index[i,2] 
        sys.cellindex[i][6,3]=cell_index_plus[i,3]   

        sys.cellindex[i][7,1]=cell_index_plus[i,1]
        sys.cellindex[i][7,2]=cell_index_plus[i,2]
        sys.cellindex[i][7,3]=cell_index[i,3] 

        sys.cellindex[i][8,1]=cell_index_plus[i,1]
        sys.cellindex[i][8,2]=cell_index_plus[i,2]
        sys.cellindex[i][8,3]=cell_index_plus[i,3]
    end
end

function grad_field!(sys::System)
    for i in 1:sys.num_atomtype
        for icell=0:sys.mesh.num_cells-2
            for jcell=0:sys.mesh.num_cells-2
                for kcell=0:sys.mesh.num_cells-2

                    icell_plus=pbc_mesh(icell+1,sys.mesh.num_cells)+1
                    jcell_plus=pbc_mesh(jcell+1,sys.mesh.num_cells)+1
                    kcell_plus=pbc_mesh(kcell+1,sys.mesh.num_cells)+1

                    icell_minus=pbc_mesh(icell-1,sys.mesh.num_cells)+1
                    jcell_minus=pbc_mesh(jcell-1,sys.mesh.num_cells)+1
                    kcell_minus=pbc_mesh(kcell-1,sys.mesh.num_cells)+1

                    icell_=icell+1
                    jcell_=jcell+1
                    kcell_=kcell+1
                    
                    g_x=sys.global_field[i][SerializeCellIndex(icell_plus,jcell_,kcell_,sys.mesh.num_cells)]-sys.global_field[i][SerializeCellIndex(icell_minus,jcell_,kcell_,sys.mesh.num_cells)]
                    g_y=sys.global_field[i][SerializeCellIndex(icell_,jcell_plus,kcell_,sys.mesh.num_cells)]-sys.global_field[i][SerializeCellIndex(icell_,jcell_minus,kcell_,sys.mesh.num_cells)]
                    g_z=sys.global_field[i][SerializeCellIndex(icell_,jcell_,kcell_plus,sys.mesh.num_cells)]-sys.global_field[i][SerializeCellIndex(icell_,jcell_,kcell_minus,sys.mesh.num_cells)]
                    
                    sys.field_gradient[i][SerializeCellGradIndex(icell_,jcell_,kcell_,1,sys.mesh.num_cells)] += 0.5*g_x
                    sys.field_gradient[i][SerializeCellGradIndex(icell_,jcell_,kcell_,2,sys.mesh.num_cells)] += 0.5*g_y
                    sys.field_gradient[i][SerializeCellGradIndex(icell_,jcell_,kcell_,3,sys.mesh.num_cells)] += 0.5*g_z
                    
                end
            end
        end
    end
end


function force_field(sys::System,field::CanonicalField,comm::MPI.Comm)
    forces=zeros(sys.local_N,3)
    if sys.first_step || sys.current_step % field.τ==0
        clear_field!(sys)
        particle2cell!(sys)
        update_cellindex!(sys)
        particle2density!(sys)
        for i in 1:sys.num_atomtype
            MPI.Allreduce!(sys.local_field[i],sys.global_field[i],sum_mesh,comm)
        end
        grad_field!(sys)

    end

    particle2cell!(sys)
    update_cellindex!(sys)
    for i in 1:sys.local_N
        sum_graident_field=zeros(3)
        for atomtype in 1:sys.num_atomtype
            field_gradient=zeros(3)
            for j in 1:8
                field_gradient[1]+=sys.field_gradient[atomtype][SerializeCellGradIndex(sys.cellindex[i][j,1],sys.cellindex[i][j,2],sys.cellindex[i][j,3],1,sys.mesh.num_cells)] * sys.particle2cell[i,j]
                field_gradient[2]+=sys.field_gradient[atomtype][SerializeCellGradIndex(sys.cellindex[i][j,1],sys.cellindex[i][j,2],sys.cellindex[i][j,3],2,sys.mesh.num_cells)] * sys.particle2cell[i,j]
                field_gradient[3]+=sys.field_gradient[atomtype][SerializeCellGradIndex(sys.cellindex[i][j,1],sys.cellindex[i][j,2],sys.cellindex[i][j,3],3,sys.mesh.num_cells)] * sys.particle2cell[i,j]
            end
            forces[i,1]+=-field.χ[sys.types[i],atomtype].*field_gradient[1]
            forces[i,2]+=-field.χ[sys.types[i],atomtype].*field_gradient[2]
            forces[i,3]+=-field.χ[sys.types[i],atomtype].*field_gradient[3]

            sum_graident_field.+=field_gradient
        end
        forces[i,1] += -1.0/field.κ * sum_graident_field[1]
        forces[i,2] += -1.0/field.κ * sum_graident_field[2]
        forces[i,3] += -1.0/field.κ * sum_graident_field[3]
    end
    
    return forces
end


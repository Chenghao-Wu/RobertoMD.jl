
function sum_mesh(point1::Float64,point2::Float64)
    return point1 + point2
end

abstract type Mesh end

struct NoMesh <:Mesh end

struct UniformMesh <:Mesh
    num_cells::Int64
    simulation_length::Float64
    cell_length::Float64
    V_cell::Float64
    mesh_index::Array{Array{Int64,1}}
    function UniformMesh(box::Array{Float64,1},Lcell::Float64)
        num_cells=convert(Int64,div(box[1],Lcell))
        real_Lcell=box[1]/num_cells
        simulation_length=box[1]
        V_cell=real_Lcell*real_Lcell*real_Lcell
        mesh_index=Array{Array{Int64,1},1}(undef,0)
        for icell=0:num_cells-2
            for jcell=0:num_cells-2
                for kcell=0:num_cells-2 
                    push!(mesh_index,[icell,jcell,kcell])
                end
            end
        end

        @info "Create UniformMesh with $(num_cells) Cells, Real Cell Length: $(real_Lcell)"
        new(num_cells::Int64,
            simulation_length::Float64,
            real_Lcell::Float64,
            V_cell::Float64,
            mesh_index::Array{Array{Int64,1}}
            )
    end
end

abstract type Field end

Base.broadcastable(x::Field) = Ref(x)

struct NoField <:Field end

struct CanonicalField <:Field

    χ::Array{Float64,2}
    κ::Float64
    τ::Int64
    Lcell::Float64
    mesh::Mesh
    local_field::Array{Array{Float64,1},1}
    global_field::Array{Array{Float64,1},1}
    field_gradient::Array{Array{Float64,1},1}
    particle2cell::Array{Array{Float64,1},1}
    cellindex::Array{Array{Int64,2},1}

    function CanonicalField(;χ::Array{Float64,2}=zeros(0,0),
                             κ::Float64=0.0,
                             τ::Int64=1,
                             Lcell::Float64=0.0,
                             mesh::Mesh=NoMesh())

        local_field=Array{Array{Float64,1},1}(undef,0)
        global_field=Array{Array{Float64,1},1}(undef,0)
        field_gradient=Array{Array{Float64,1},1}(undef,0)
        
        particle2cell=Array{Array{Float64,1},1}(undef,0)
        cellindex=Array{Array{Int64,2},1}(undef,0)

        new(χ::Array{Float64,2},
            κ::Float64,
            τ::Int64,
            Lcell::Float64,
            mesh::Mesh,
            local_field::Array{Array{Float64,1},1},
            global_field::Array{Array{Float64,1},1},
            field_gradient::Array{Array{Float64,1},1},
            particle2cell::Array{Array{Float64,1},1},
            cellindex::Array{Array{Int64,2},1})
    end
end



function init(::CanonicalField,input::Dict,configuration::Configuration)
    num_atomtypes=length(keys(input["Canonical field"]["χ"]))
    chi=zeros(num_atomtypes,num_atomtypes)
    for atomi in 1:num_atomtypes
        chi_i=convert(Array{Float64,1},input["Canonical field"]["χ"][string(atomi)])
        for atomj in 1:num_atomtypes
            chi[atomi,atomj]=chi_i[atomj]
        end
    end

    if "uniform mesh" in keys(input["Canonical field"])
        mesh=UniformMesh(configuration.box,input["Canonical field"]["Lcell"])
    end

    field=CanonicalField(χ=chi,
                         κ=input["Canonical field"]["κ"],
                         τ=input["Canonical field"]["update"],
                         Lcell=input["Canonical field"]["Lcell"],
                         mesh=mesh)

    for i in 1:num_atomtypes
        push!(field.local_field,zeros(mesh.num_cells*mesh.num_cells*mesh.num_cells))
        push!(field.global_field,zeros(mesh.num_cells*mesh.num_cells*mesh.num_cells))
        push!(field.field_gradient,zeros(mesh.num_cells*mesh.num_cells*mesh.num_cells*3))
    end
    return field
end

function sum_mesh(point1::Float64,point2::Float64)
    return point1 + point2
end

abstract type Field end

struct NoField <:Field end

struct CanonicalField <:Field
    χ::Array{Float64,2}
    κ::Float64
    τ::Int64
    Lcell::Float64
    function CanonicalField(χ::Array{Float64,2},
                            κ::Float64,
                            τ::Int64,
                            Lcell::Float64)
        new(χ::Array{Float64,2},
            κ::Float64,
            τ::Int64,
            Lcell::Float64)
    end
end

abstract type Mesh end

struct NoMesh <:Mesh end

struct UniformMesh <:Mesh
    num_cells::Int64
    simulation_length::Float64
    cell_length::Float64
    function UniformMesh(pbc::PBC,Lcell::Float64)
        num_cells=convert(Int64,div(pbc.Lx,Lcell))
        real_Lcell=pbc.Lx/num_cells
        simulation_length=pbc.Lx
        @info "Create UniformMesh with $(num_cells) Cells, Real Cell Length: $(real_Lcell)"
        new(num_cells::Int64,
            simulation_length::Float64,
            real_Lcell::Float64
            )
    end
end


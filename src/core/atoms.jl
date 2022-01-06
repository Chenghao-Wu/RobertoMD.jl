
struct Atoms

    total_N::Array{Int64,1}
    local_N::Array{Int64,1}

    rho0::Array{Float64,1}

    num_atomtype::Array{Int64,1}

    types::Array{Int64,1}
    types_all::Array{Int64,1}

    masses::Array{Float64,1}
    masses_all::Array{Float64,1}

    coords::Vector{Array{Float64,1}}
    coords_all::Vector{Array{Float64,1}}

    vels::Vector{Array{Float64,1}}
    vels_all::Vector{Array{Float64,1}}

    forces::Vector{Array{Float64,1}}
    forces_all::Vector{Array{Float64,1}}

    energies::Array{Float64,1}
    energies_all::Array{Float64,1}

    bonds::Vector{Array{Int64,1}}
    bonds_all::Vector{Array{Int64,1}}

    function Atoms()
        total_N=[0]
        local_N=[0]
        rho0=[0.0]
        num_atomtype=[0]

        types=Vector{Int64}(undef,0)
        types_all=Vector{Int64}(undef,0)

        masses=Vector{Float64}(undef,0)
        masses_all=Vector{Float64}(undef,0)

        coords=Array{Array{Float64,1},1}(undef,0)
        coords_all=Array{Array{Float64,1},1}(undef,0)

        vels=Array{Array{Float64,1},1}(undef,0)
        vels_all=Array{Array{Float64,1},1}(undef,0)

        forces=Array{Array{Float64,1},1}(undef,0)
        forces_all=Array{Array{Float64,1},1}(undef,0)

        energies=Vector{Float64}(undef,0)
        energies_all=Vector{Float64}(undef,0)

        bonds=Vector{Vector{Int64}}(undef,0)
        bonds_all=Vector{Vector{Int64}}(undef,0)

        new(total_N::Array{Int64,1},
            local_N::Array{Int64,1},
            rho0::Array{Float64,1},
            num_atomtype::Array{Int64,1},
            types::Array{Int64,1},
            types_all::Array{Int64,1},
        
            masses::Array{Float64,1},
            masses_all::Array{Float64,1},
        
            coords::Vector{Array{Float64,1}},
            coords_all::Vector{Array{Float64,1}},
        
            vels::Vector{Array{Float64,1}},
            vels_all::Vector{Array{Float64,1}},
        
            forces::Vector{Array{Float64,1}},
            forces_all::Vector{Array{Float64,1}},
        
            energies::Array{Float64,1},
            energies_all::Array{Float64,1},
        
            bonds::Vector{Array{Int64,1}},
            bonds_all::Vector{Array{Int64,1}})
    end
end

function init_atoms!(atoms::Atoms,
                    types_all::Array{Int64,1},
                    masses_all::Array{Float64,1},
                    coords_all::Array{Float64,2},
                    forces_all::Array{Float64,2},
                    vels_all::Array{Float64,2},
                    energies_all::Array{Float64,1})
    for i in 1:atoms.total_N[1]
        push!(atoms.types_all,types_all[i])
        push!(atoms.masses_all,masses_all[i])
        push!(atoms.coords_all,coords_all[i,1:3])
        push!(atoms.forces_all,forces_all[i,1:3])
        push!(atoms.vels_all,vels_all[i,1:3])
        push!(atoms.energies_all,energies_all[i])
    end
end


function init_bonds!(atoms::Atoms,
                    bonds_all::Array{Int64,2})
    for bondi in 1:size(bonds_all)[1]
        push!(atoms.bonds_all,bonds_all[bondi,1:3])
    end
end
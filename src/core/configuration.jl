export Configuration,ReadConfig

function ReadConfig(filename::String)
    input = open(filename, "r")
    inputs=JSON.parse(input)
    close(input)
    return inputs
end

mutable struct Molecule
    type::Array{Int64,1}
    mass::Array{Float64,1}
    coords::Array{Float64,2}
    bonds::Array{Int64,2}
    angles::Array{Int64,2}
    vels::Array{Float64,2}
end

#function Molecule(type::Array{Int64,1},
#    mass::Array{Float64,1},
#    coords::Array{Float64,2},
#    bonds::Array{Array{Int64,2},1},
#    angles::Array{Array{Int64,2},1},
#    vels::Array{Float64,2})
#new(
#type::Array{Int64,1},
#mass::Array{Float64,1},
#coords::Array{Float64,2},
#bonds::Array{Array{Int64,2},1},
#angles::Array{Array{Int64,2},1},
#vels::Array{Float64,2}
#)
#end

struct Configuration
    start::Int64
    box::Array{Float64,1}
    moles::Array{Molecule,1}
    moles_length::Array{Int64,1}
    function Configuration(data::Dict)
        start=1
        box=zeros(3)
        if "box" in keys(data)
            box=convert(Array{Float64,1}, data["box"])
            @info "Box Lx Ly Lz: $(box[1]) $(box[2]) $(box[3])"
        end

        moles=Vector{Molecule}(undef,0)
        moles_length=zeros(Int64,0)

        if "molecules" in keys(data)
            molecules=data["molecules"]
            num_moles=length(keys(molecules))

            moles=Vector{Molecule}(undef,num_moles)

            moles_length=zeros(Int64,num_moles)
            for molei in 1:num_moles
                moles_length[molei]=length(keys(molecules[string(molei)]["atoms"]))
                #molecule=Molecule(moles_length[molei])
                type=zeros(Int64,moles_length[molei])
                mass=zeros(Int64,moles_length[molei])
                bonds=zeros(Int64,0,3)
                angles=zeros(Int64,moles_length[molei],4)
                coords=zeros(Float64,moles_length[molei],3)
                vels=zeros(Float64,moles_length[molei],3)
                for atomi in 1:moles_length[molei]
                    type[atomi]=molecules[string(molei)]["atoms"][string(atomi)]["type"]
                    mass[atomi]=molecules[string(molei)]["atoms"][string(atomi)]["mass"]
                    coords[atomi,1:3]=convert(Array{Float64,1},molecules[string(molei)]["atoms"][string(atomi)]["coords"])
                    if "vels" in keys(molecules[string(molei)]["atoms"][string(atomi)])
                        vels[atomi,1:3]=convert(Array{Float64,1},molecules[string(molei)]["atoms"][string(atomi)]["vels"])
                    end
                end
                if "bonds" in keys(molecules[string(molei)])
                    num_bonds=length(keys(molecules[string(molei)]["bonds"]))
                    bonds=zeros(Int64,num_bonds,3)
                    for bondi in 1:num_bonds
                        bonds[bondi,1:3]=convert(Array{Int64,1},molecules[string(molei)]["bonds"][string(bondi)])
                    end
                end
                moles[molei]=Molecule(type,mass,coords,bonds,angles,vels)
            end
        end
        
        new(start,box,
            moles,
            moles_length)
    end
end
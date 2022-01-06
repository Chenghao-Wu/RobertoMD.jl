export Harmonic

abstract type BdInter end
Base.broadcastable(x::BdInter) = Ref(x)
struct NoBond<:BdInter end

struct Harmonic<:BdInter
    k::Array{Float64,1}
    l0::Array{Float64,1}
    function Harmonic(;k::Array{Float64,1}=[0.0],
                        l0::Array{Float64,1}=[0.0])
        new(k::Array{Float64,1},
            l0::Array{Float64,1})
    end
end

function init(::Harmonic,input::Dict)
    k_matrix=zeros(length(keys(input["harmonic bond"]["k"])))
    l0_matrix=zeros(length(keys(input["harmonic bond"]["l0"])))

    for atomi in 1:length(keys(input["harmonic bond"]["k"]))
        k_matrix[atomi]=input["harmonic bond"]["k"][string(atomi)]
        l0_matrix[atomi]=input["harmonic bond"]["l0"][string(atomi)]
    end
    bdinter=Harmonic(k=k_matrix,l0=l0_matrix)
    return bdinter
end
export Harmonic

abstract type BdInter end
Base.broadcastable(x::BdInter) = Ref(x)
struct NoBond<:BdInter end

struct Harmonic<:BdInter
    k::Array{Float64,1}
    l0::Array{Float64,1}
    function Harmonic(  k::Array{Float64,1},
                        l0::Array{Float64,1})
        new(k::Array{Float64,1},
            l0::Array{Float64,1})
    end
end
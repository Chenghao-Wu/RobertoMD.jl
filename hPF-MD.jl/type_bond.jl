export Harmonic

abstract type BdInter end

struct NoBond<:BdInter end

struct Harmonic<:BdInter
    k::Float64
    l0::Float64
    function Harmonic(k::Float64,l0::Float64)
        new(k::Float64,l0::Float64)
    end
end
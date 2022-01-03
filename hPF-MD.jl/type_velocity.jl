abstract type Velocity end

struct NoVelocity<:Velocity end

struct GaussianVelocity<:Velocity
    T::Float64
    seed::Int64
    function GaussianVelocity(  T::Float64,
                                seed::Int64)
        new(T::Float64,
            seed::Int64)
    end
end

struct ZeroVelocity<:Velocity end
export LangevinNVT

abstract type Thermostat end

struct NoThermostat<:Thermostat end

struct LangevinNVT<:Thermostat
    gamma::Float64
    invdt::Float64
    function LangevinNVT(gamma::Float64,dt::Float64)
        invdt=1.0/dt
        new(gamma::Float64,invdt::Float64)
    end
end

struct BerendsenNVT <: Thermostat 
    τ::Float64
end
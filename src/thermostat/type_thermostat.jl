export LangevinNVT

abstract type Thermostat end

struct NoThermostat<:Thermostat end

struct LangevinNVT<:Thermostat
    invdt::Float64
    gamma::Float64
    function LangevinNVT(;dt::Float64=0.0,gamma::Float64=0.0)
        invdt=1.0/dt
        new(invdt::Float64,
            gamma::Float64)
    end
end

function init(::LangevinNVT,input::Dict)
    gamma=input["LangevinNVT"]["gamma"]
    dt=input["dt"]
    thermostat=LangevinNVT(dt=dt,gamma=gamma)
    return thermostat
end


struct BerendsenNVT <: Thermostat 
    τ::Float64
    function BerendsenNVT(;τ=0.0)
        new(τ::Float64)
    end
end

function init(::BerendsenNVT,input::Dict)
    tau=input["BerendsenNVT"]["tau"]
    thermostat=BerendsenNVT(τ=tau)
    return thermostat
end
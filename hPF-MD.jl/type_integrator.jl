abstract type Integrator end

struct NoIntegrator <: Integrator end

struct VVIntegrator <: Integrator 
    hdtsqr::Float64
    function VVIntegrator(dt::Float64)
        hdtsqr=dt*dt/2.0 
        new(hdtsqr::Float64)
    end
end # Velocity Verlet
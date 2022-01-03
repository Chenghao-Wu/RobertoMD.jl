function UpdatePosition!(sys::System,integrator::NoIntegrator)
end

function UpdateVelocity!(sys::System,integrator::NoIntegrator)
end

function UpdatePosition!(sys::System,integrator::VVIntegrator)
    sys.coords .+= sys.dt.*sys.vels.+integrator.hdtsqr .* sys.forces ./ sys.masses
end

function UpdateVelocity!(sys::System,integrator::VVIntegrator,comm::MPI.Comm)
    sys.vels .+= 0.5 .* sys.dt .* sys.forces ./ sys.masses
end
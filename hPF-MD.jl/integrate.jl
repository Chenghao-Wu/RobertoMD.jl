function UpdatePosition!(sys::System,integrator::NoIntegrator)
end

function UpdateVelocity!(sys::System,integrator::NoIntegrator)
end

function UpdatePosition!(sys::System,integrator::VVIntegrator)
    #sys.coords .+= sys.dt.*sys.vels.+integrator.hdtsqr .* sys.forces ./ sys.masses
    UpdatePosition_particle!.(sys.coords,sys.vels,sys.forces,sys.masses,sys.dt,integrator)
end

function UpdatePosition_particle!(  coords::Array{Float64,1},
                                    vels::Array{Float64,1},
                                    forces::Array{Float64,1},
                                    masses::Float64,
                                    dt::Float64,
                                    integrator::VVIntegrator)
    coords[1] += dt* vels[1] + integrator.hdtsqr * forces[1] / masses
    coords[2] += dt* vels[2] + integrator.hdtsqr * forces[2] / masses
    coords[3] += dt* vels[3] + integrator.hdtsqr * forces[3] / masses
end

function UpdateVelocity!(sys::System,integrator::VVIntegrator,comm::MPI.Comm)
    #sys.vels .+= 0.5 .* sys.dt .* sys.forces ./ sys.masses
    UpdateVelocity_particle!.(sys.vels,sys.forces,sys.masses,sys.dt,integrator)
end

function UpdateVelocity_particle!(  vels::Array{Float64,1},
                                    forces::Array{Float64,1},
                                    masses::Float64,
                                    dt::Float64,
                                    integrator::VVIntegrator)
    vels[1] += 0.5 * dt * forces[1] / masses
    vels[2] += 0.5 * dt * forces[2] / masses
    vels[3] += 0.5 * dt * forces[3] / masses
end
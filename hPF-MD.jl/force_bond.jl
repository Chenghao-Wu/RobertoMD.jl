function force_bond!(sys::System,bdinter::NoBond)
    forces=zeros(sys.local_N,3)
    return forces
end

function force_bond!(sys::System,bdinter::Harmonic)
    force_bond_atomic!.(Ref(sys.forces),sys.bonds,Ref(sys.coords),bdinter,sys.pbc.Lx,sys.pbc.Ly,sys.pbc.Lz)
end

function force_bond_atomic!(forces::Array{Array{Float64,1}},
                            bond::Array{Int64,1},
                            coords::Array{Array{Float64,1}},
                            bdinter::Harmonic,lx::Float64,ly::Float64,lz::Float64)
    
    delta_x = coords[bond[2]][1]- coords[bond[3]][1]
    delta_y = coords[bond[2]][2]- coords[bond[3]][2]
    delta_z = coords[bond[2]][3]- coords[bond[3]][3]

    # consider pbc
    delta_x -= lx * RInt( delta_x / lx)
    delta_y -= ly * RInt( delta_y / ly)
    delta_z -= lz * RInt( delta_z / lz)

    r=(delta_x^2+delta_y^2+delta_z^2)^0.5

    κ=bdinter.k[bond[1]] #κ
    l0=bdinter.l0[bond[1]] #r_0

    force_=κ*(l0/r-1.0)
    force_divr = force_/r

    forces[bond[2]][1] += delta_x * force_divr
    forces[bond[2]][2] += delta_y * force_divr
    forces[bond[2]][3] += delta_z * force_divr

    forces[bond[3]][1] += -delta_x * force_divr
    forces[bond[3]][2] += -delta_y * force_divr
    forces[bond[3]][3] += -delta_z * force_divr
    #@show bond[2],bond[3],force_divr,delta_x,delta_y,delta_z
end
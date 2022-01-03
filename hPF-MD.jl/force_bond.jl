function force_bond(sys::System,bdinter::NoBond)
    forces=zeros(sys.local_N,3)
    return forces
end

function force_bond(sys::System,bdinter::Harmonic)
    forces=zeros(sys.local_N,3)
    return forces
end

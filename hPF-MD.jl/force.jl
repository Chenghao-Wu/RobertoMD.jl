

function CalculateForce!(sys::System,comm::MPI.Comm)
    if sys.bdinter!=NoBond()
        sys.forces.+=force_bond(sys,sys.bdinter)
    end
    if sys.field!=NoField()
        sys.forces.+=force_field(sys,sys.field,comm)
    end
end

function ClearForce!(sys::System)
    fill!(sys.forces,0.0)
end
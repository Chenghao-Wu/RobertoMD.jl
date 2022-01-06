

function CalculateForce!(sys::System,comm::MPI.Comm)
    for bond in sys.bdinter
        if bond!=NoBond()
            force_bond!(sys,bond)
        end
    end

    if sys.field!=NoField()
        force_field!(sys,sys.field,comm)
    end
    
end

function ClearForce!(sys::System)
    fill!.(sys.atoms.forces,0.0)
end
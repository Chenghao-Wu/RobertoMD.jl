
function DeSerializeAllCoord!(sys::System,array::Array{Float64,1})
    for i in 0:div(size(array)[1],3)-1
        sys.atoms.coords_all[i+1][1]=array[i*3+1]
        sys.atoms.coords_all[i+1][2]=array[i*3+2]
        sys.atoms.coords_all[i+1][3]=array[i*3+3]
    end
end

function DeSerializeAllVels!(sys::System,array::Array{Float64,1})
    for i in 0:div(size(array)[1],3)-1
        sys.atoms.vels_all[i+1][1]=array[i*3+1]
        sys.atoms.vels_all[i+1][2]=array[i*3+2]
        sys.atoms.vels_all[i+1][3]=array[i*3+3]
    end
end

function DeSerializeAllBonds!(sys::System,array::Array{Int64,1})
    for i in 0:div(size(array)[1],3)-1
        sys.atoms.bonds_all[i+1][1]=array[i*3+1]
        sys.atoms.bonds_all[i+1][2]=array[i*3+2]
        sys.atoms.bonds_all[i+1][3]=array[i*3+3]
    end
end

function DeSerializeAllMasses!(sys::System,array::Array{Float64,1})
    for i in 0:div(size(array)[1],1)-1
        sys.atoms.masses_all[i+1]=array[i+1]
    end
end

function DeSerializeAllTypes!(sys::System,array::Array{Int64,1})
    for i in 0:div(size(array)[1],1)-1
        sys.atoms.types_all[i+1]=array[i+1]
    end
end

function DeInitMolecules(sys::System)
    molecules=Dict()
    index_atom=1
    index_bond=1
    for comm_mole_i in sys.mpi.mole_index_comms
        for mole_i in comm_mole_i

            molecules[string(mole_i)]=Dict()

            molecules[string(mole_i)]["bonds"]=Dict()
            bond_mole_length=sys.mpi.bond_moles_length[mole_i]
            
            for bondi in 1:bond_mole_length
                bond_=zeros(Int64,3)
                bond_[1]=sys.atoms.bonds_all[index_bond][1]
                bond_[2]=sys.atoms.bonds_all[index_bond][2] - index_atom + 1
                bond_[3]=sys.atoms.bonds_all[index_bond][3] - index_atom + 1
                molecules[string(mole_i)]["bonds"][string(bondi)]=bond_
                index_bond=index_bond+1
            end

            molecules[string(mole_i)]["atoms"]=Dict()
            mole_length=sys.mpi.moles_length[mole_i]
            for atomi in 1:mole_length
                molecules[string(mole_i)]["atoms"][string(atomi)]=Dict( "type"=>sys.atoms.types_all[index_atom],
                                                                        "mass"=>sys.atoms.masses_all[index_atom],
                                                                        "coords"=>sys.atoms.coords_all[index_atom],
                                                                        "vels"=>sys.atoms.vels_all[index_atom])
                index_atom=index_atom+1
            end
            
        end
    end
    return molecules
end

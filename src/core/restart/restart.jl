function DeInitInformation!(sys::System,comm::MPI.Comm,root::Int64)
    coords_all=SerializeArray(sys.atoms.coords_all)
    MPI.Gatherv!(SerializeArray(sys.atoms.coords),VBuffer(coords_all,sys.mpi.commsizes*3),root,comm)
    DeSerializeAllCoord!(sys,coords_all)

    vels_all=SerializeArray(sys.atoms.vels_all)
    MPI.Gatherv!(SerializeArray(sys.atoms.vels),VBuffer(vels_all,sys.mpi.commsizes*3),root,comm)
    DeSerializeAllVels!(sys,vels_all)

    bonds_all=SerializeArray(sys.atoms.bonds_all)
    MPI.Gatherv!(SerializeArray(sys.atoms.bonds),VBuffer(bonds_all,sys.mpi.bond_commsizes*3),root,comm)
    DeSerializeAllBonds!(sys,bonds_all)

    types_all=SerializeArray(sys.atoms.types_all)
    MPI.Gatherv!(SerializeArray(sys.atoms.types),VBuffer(types_all,sys.mpi.commsizes),root,comm)
    DeSerializeAllTypes!(sys,types_all)

    masses_all=SerializeArray(sys.atoms.masses_all)
    MPI.Gatherv!(SerializeArray(sys.atoms.masses),VBuffer(masses_all,sys.mpi.commsizes),root,comm)
    DeSerializeAllMasses!(sys,masses_all)
end

function RestartInformation(sys::System,comm::MPI.Comm,root::Int64)
    MPI.Barrier(comm)
    if sys.restart!=NoRestart()
        if sys.restart.freq==-1 && sys.current_step[1]==sys.steps+sys.start_step-1
            DeInitInformation!(sys,comm,root)
            if MPI.Comm_rank(comm)==root
                Restart(sys,sys.restart)
            end
        elseif sys.current_step[1]%sys.restart.freq==0 || sys.current_step[1]==sys.steps+sys.start_step-1
            DeInitInformation!(sys,comm,root)
            if MPI.Comm_rank(comm)==root
                Restart(sys,sys.restart)
            end
        end
    end
    MPI.Barrier(comm)
end

function Restart(sys::System,restart::JSONRestart)

    molecules=DeInitMolecules(sys)

    configs=Dict(
        "box"=>[sys.pbc.Lx,sys.pbc.Ly,sys.pbc.Lz],
        "molecules"=>molecules
        )
    #control_file=open("control.json","w")
    JSON.print(restart.io,configs,4)
end
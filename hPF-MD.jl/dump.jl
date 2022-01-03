
function DumpInformation(sys::System,comm::MPI.Comm,root::Int64)
    MPI.Barrier(comm)
    
    if sys.trjdump!=NoDump()
        coords_all=SerializeArray(sys.coords_all)
        MPI.Gatherv!(SerializeArray(sys.coords),VBuffer(coords_all,sys.commsizes*3),root,comm)
        sys.coords_all=DeSerializeArray(coords_all)
        if MPI.Comm_rank(comm)==root
            DumpTrajectory(sys,sys.trjdump)
        end
    end
    MPI.Barrier(comm)
end

function DumpTrajectory(sys::System,dump::LAMMPSTrj)
    if sys.first_step || sys.current_step%dump.freq==0
        string_out="ITEM: TIMESTEP"
        string_out*="\n"
        string_out*=string(sys.current_step)
        string_out*="\n"
        string_out*="ITEM: NUMBER OF ATOMS"
        string_out*="\n"
        string_out*=string(sys.total_N)
        string_out*="\n"
        string_out*="ITEM: BOX BOUNDS pp pp pp"
        string_out*="\n"
        string_out*="0 "*string(sys.pbc.Lx)
        string_out*="\n"
        string_out*="0 "*string(sys.pbc.Ly)
        string_out*="\n"
        string_out*="0 "*string(sys.pbc.Lz)
        string_out*="\n"

        string_out*="ITEM: ATOMS"
        string_out*=" id"
        string_out*=" type x y z"

        if dump.is_velocity==true
            string_out*=" vx vy vz"
        end
        if dump.is_force==true
            string_out*=" fx fy fz"
        end
        string_out*="\n"
        write(dump.io,string_out)
        for atomii=1:sys.total_N
            @fastmath @inbounds begin
            string_out=""
            string_out*=string(atomii)*" "
            string_out*=string(sys.types_all[atomii])*" "*string(sys.coords_all[atomii,1])*" "*string(sys.coords_all[atomii,2])*" "*string(sys.coords_all[atomii,3])*" "

            if dump.is_velocity==true
                string_out*=string(sys.vels[atomii])*" "*string(sys.vels[atomii])*" "*string(sys.vels[atomii])*" "
            end
            if dump.is_force==true
                string_out*=string(sys.forces[atomii])*" "*string(sys.forces[atomii])*" "*string(sys.forces[atomii])*" "
            end
            string_out*="\n"
            write(dump.io,string_out)
            end
        end
    end
end
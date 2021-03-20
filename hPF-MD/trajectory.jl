

export XYZTrajectoryDump,LAMMPSTrajectoryDump

function apply_dump(fileO,step::Int64,args::system,atoms::Vector{atom},forces::Vector{force},velocities::Vector{velocity},::NoTrajectoryDump) end

struct LAMMPSTrajectoryDump <: TrajectoryDump
    velocity::Bool
    force::Bool
    image::Bool
    id::Bool
end

function apply_dump(fileO,step::Int64,args::system,atoms::Vector{atom},forces::Vector{force},velocities::Vector{velocity},dump::LAMMPSTrajectoryDump) 
    string_out="ITEM: TIMESTEP"
    string_out*="\n"
    string_out*=string(step)
    string_out*="\n"
    string_out*="ITEM: NUMBER OF ATOMS"
    string_out*="\n"
    string_out*=string(length(atoms))
    string_out*="\n"
    string_out*="ITEM: BOX BOUNDS pp pp pp"
    string_out*="\n"
    string_out*="0 "*string(args.box[1])
    string_out*="\n"
    string_out*="0 "*string(args.box[2])
    string_out*="\n"
    string_out*="0 "*string(args.box[3])
    string_out*="\n"

    string_out*="ITEM: ATOMS"
    if dump.id==true
        string_out*=" id"
    end
    string_out*=" type x y z"
    if dump.image==true
        string_out*=" ix iy iz"
    end
    if dump.velocity==true
        string_out*=" vx vy vz"
    end
    if dump.force==true
        string_out*=" fx fy fz"
    end
    string_out*="\n"
    write(fileO,string_out)

    for atomii=1:length(atoms)
        string_out=""
        if dump.id==true
            string_out*=string(atomii)*" "
        end
        string_out*=string(atoms[atomii].type)*" "*string(atoms[atomii].pos[1])*" "*string(atoms[atomii].pos[2])*" "*string(atoms[atomii].pos[3])*" "
        if dump.image==true
            string_out*=string(atoms[atomii].image[1])*" "*string(atoms[atomii].image[2])*" "*string(atoms[atomii].image[3])*" "
        end
        if dump.velocity==true
            string_out*=string(velocities[atomii].velocity[1])*" "*string(velocities[atomii].velocity[2])*" "*string(velocities[atomii].velocity[3])*" "
        end
        if dump.force==true
            string_out*=string(forces[atomii].force[1])*" "*string(forces[atomii].force[2])*" "*string(forces[atomii].force[3])*" "
        end
        string_out*="\n"
        write(fileO,string_out)
    end

end

struct XYZTrajectoryDump <: TrajectoryDump
end

function apply_dump(fileO,step::Int64,args::system,atoms::Vector{atom},forces::Vector{force},velocities::Vector{velocity},::XYZTrajectoryDump) 
    string_out=""
    string_out=string(length(atoms))
    string_out*="\n"
    write(fileO,string_out)
    string_out=""
    string_out*=string(step)
    string_out*="\n"
    write(fileO,string_out)
    for atomii=1:length(atoms)
        string_out=""
        type_=atoms[atomii].type
        string_out=string(type_)*" "*string(atoms[atomii].pos[1])*" "*string(atoms[atomii].pos[2])*" "*string(atoms[atomii].pos[3])
        string_out*="\n"
        write(fileO,string_out)
    end
end
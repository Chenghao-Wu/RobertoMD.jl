

export XYZTrajectoryDump

function apply_dump(atoms::Vector{atom},forces::Vector{force},velocities::Vector{velocity},::NoTrajectoryDump) end

struct LAMMPSTrajectoryDump <: TrajectoryDump
    velocity::Bool
    force::Bool
end

function apply_dump(fileO,atoms::Vector{atom},forces::Vector{force},velocities::Vector{velocity},::LAMMPSTrajectoryDump) 
end

struct XYZTrajectoryDump <: TrajectoryDump
end

function apply_dump(fileO,step::Int64,atoms::Vector{atom},forces::Vector{force},velocities::Vector{velocity},::XYZTrajectoryDump) 
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


abstract type Dump end

struct NoDump<: Dump end

struct LAMMPSTrj <: Dump
    filename::String
    freq::Int64
    is_velocity::Bool
    is_force::Bool
    io::IOStream
    function LAMMPSTrj( filename::String,
                        freq::Int64;
                        is_velocity=false,
                        is_force=false)
        io=open(filename,"w")
        new(filename::String,
            freq::Int64,
            is_velocity::Bool,
            is_force::Bool,
            io::IOStream
        )
    end
end
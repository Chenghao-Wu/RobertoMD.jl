

abstract type Restart end

struct NoRestart<: Restart end

struct JSONRestart <: Restart
    filename::String
    freq::Int64
    io::IOStream
    function JSONRestart(filename::String,
                        freq::Int64)
        io=open(filename,"w")
        new(filename::String,
            freq::Int64,
            io::IOStream
        )
    end
end
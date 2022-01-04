abstract type Logger end

struct NoLogger <:Logger end

struct ThermoLogger <:Logger
    Step::Bool
    Temp::Bool
    Energy::Bool
    Momentum::Bool
    index::Array{Int64,1}
    Write::Bool
    File::String
    function ThermoLogger(; Energy=false,
                            Momentum=false,
                            Write=false,
                            File="" )
                            Step=true
                            Temp=true
                            index=[1]
                            new(
                                Step::Bool,
                                Temp::Bool,
                                Energy::Bool,
                                Momentum::Bool,
                                index::Array{Int64,1},
                                Write::Bool,
                                File::String
                            )
    end
end


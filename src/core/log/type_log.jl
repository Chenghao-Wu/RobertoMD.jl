abstract type Logger end

struct NoLogger <:Logger end


struct ThermoLogger <:Logger
    freq::Int64
    Step::Bool
    Temp::Bool
    Energy::Bool
    Momentum::Bool
    index::Array{Int64,1}
    Write::Bool
    File::String
    Info::Array{Array{Float64,1}}
    function ThermoLogger(freq::Int64,
                            ; Energy=false,
                            Momentum=false,
                            Write=false,
                            File="" )
                            Step=true
                            Temp=true
                            index=[1]
                            info=Array{Array{Float64,1},1}(undef,0)
                            new(freq::Int64,
                                Step::Bool,
                                Temp::Bool,
                                Energy::Bool,
                                Momentum::Bool,
                                index::Array{Int64,1},
                                Write::Bool,
                                File::String,
                                info::Array{Array{Float64,1}}
                            )
    end
end

function init(logger::Logger,input::Dict,start_step::Int64)
    #init!()
    thermofreq=input["thermo information"]["freq"]

    @info "Log Thermo Information Every: $(thermofreq)"

    L=2
    energy_=false
    momentum_=false
    write_=false
    file_=""
    for i in 1:length(keys(input["thermo information"]))
        if "energy" in keys(input["thermo information"])
            if input["thermo information"]["energy"]
                energy_=true
                L+=1
            end
        end
        if "momentum" in keys(input["thermo information"])
            if input["thermo information"]["momentum"]
                momentum_=true
                L+=1
            end
        end
        if "write" in keys(input["thermo information"])
            write_=input["thermo information"]["write"]
        end
        if "file" in keys(input["thermo information"])
            file_=input["thermo information"]["file"]
        end
    end
    
    thermologger=ThermoLogger(thermofreq,Energy=energy_,Momentum=momentum_,Write=write_,File=file_)

    if thermofreq!=1
        push!(thermologger.Info,zeros(L))
    end
    for stepi in 1:input["steps"]
        
        if stepi%thermofreq==0.0
            push!(thermologger.Info,zeros(L))
        end
    end
    return thermologger
end
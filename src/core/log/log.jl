function LoggerInformation!(sys::System,comm::MPI.Comm,root::Int64)
    if sys.thermologger!=NoLogger()
        LoggerThermo!(sys,comm,root)
    end
end


function LoggerThermo!(sys::System,comm::MPI.Comm,root::Int64)
    if sys.first_step[1] || sys.current_step[1]%sys.thermologger.freq==0.0
        sys.thermologger.Info[sys.thermologger.index[1]][1]=sys.current_step[1]
        sys.thermologger.Info[sys.thermologger.index[1]][2]=sys.current_temp[1]
        index=3
        if sys.thermologger.Momentum
            sys.thermologger.Info[sys.thermologger.index[1]][index]=sys.current_momentum[1]
            index+=1
        end
        if MPI.Comm_rank(comm)==root
            step=sys.current_step[1]
            temp=sys.current_temp[1]
            @show step,temp
        end
        sys.thermologger.index[1]+=1
    end
    if sys.current_step[1]==sys.steps+sys.start_step-1
        if MPI.Comm_rank(comm)==root
            if sys.thermologger.Write
                write_log(sys,sys.thermologger)
            end
        end
    end
end

function write_log(sys::System,logger::ThermoLogger)
    io=open(logger.File,"w")
    string_out="Step Temp"
    if logger.Momentum
        string_out*=" Momentum"
    end
    string_out*="\n"
    write(io,string_out)
    
    for i in sys.thermologger.Info
        index=3
        string_out=""
        string_out*=string(Int(i[1]))*" "
        string_out*=string(i[2])*" "
        if sys.thermologger.Momentum
            string_out*=string(i[index])*" "
            index+=1
        end
        string_out*="\n"
        write(io,string_out)
    end
    close(io)
end
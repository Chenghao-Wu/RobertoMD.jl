
export Logger,logger_init

mutable struct Logger_ <: Logger 
    Timestep::Bool
    Temp::Bool
    PotentialEnergy::Bool
    KineticEnergy::Bool
    Momentum::Bool
end

function logger_init()
    Timestep=true
    Temp=true
    PotentialEnergy=false
    KineticEnergy=false
    Momentum=false
    return Logger_(Timestep,Temp,PotentialEnergy,KineticEnergy,Momentum)
end

function apply_logger!(fileO,firststep::Bool,step::Int64,velocities::Vector{velocity},energy::Vector{Float64},atoms::Vector{atom},nl::NoLogger) end

function apply_logger!(fileO,firststep::Bool,step::Int64,velocities::Vector{velocity},energy::Vector{Float64},atoms::Vector{atom},log::Logger_) 
    totalkineticenergy=0
    if firststep==true
        string_out=""
        if log.Timestep==true
            string_out*="timestep"*" "
        end
        if log.Temp==true
            string_out*="Temp"*" "
        end
        if log.PotentialEnergy==true
            string_out*="PotentialEnergy"*" "
        end
        if log.KineticEnergy==true
            string_out*="KineticEnergy"*" "
        end
        if log.Momentum==true
            string_out*="Momentum"
        end
        string_out*="\n"
        write(fileO,string_out)
        flush(fileO)
        string_out=""
        if log.Timestep==true
            string_out*=string(step)*" "
        end
        if log.Temp==true
            totalkineticenergy=calc_totalkineticenergy(atoms,velocities)
            inst_temp=calc_inst_temp(atoms,totalkineticenergy)
            string_out*=format(inst_temp,precision=6)*" "
        end
        if log.PotentialEnergy==true
            total_energy=sum(energy)
            string_out*=format(total_energy,precision=6)*" "
        end
        if log.KineticEnergy==true
            string_out*=format(totalkineticenergy,precision=6)*" "
        end
        if log.Momentum==true
            Momentum=calc_totalmomentum(atoms,velocities)
            string_out*=format(Momentum,precision=6)*" "
        end

        string_out*="\n"
        write(fileO,string_out)
        flush(fileO)
    else
        string_out=""
        if log.Timestep==true
            string_out*=string(step)*" "
        end
        if log.Temp==true
            totalkineticenergy=calc_totalkineticenergy(atoms,velocities)
            inst_temp=calc_inst_temp(atoms,totalkineticenergy)
            string_out*=format(inst_temp,precision=6)*" "
        end
        if log.PotentialEnergy==true
            total_energy=sum(energy)
            string_out*=format(total_energy,precision=6)*" "
        end
        if log.KineticEnergy==true
            string_out*=format(totalkineticenergy,precision=6)*" "
        end
        if log.Momentum==true
            Momentum=calc_totalmomentum(atoms,velocities)
            string_out*=format(Momentum,precision=6)*" "
        end
        string_out*="\n"
        write(fileO,string_out)
        flush(fileO)
    end
end
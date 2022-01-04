
include("/home/zwq2834/development/hPF-MD.jl/hPF-MD.jl/hPF-MD.jl")

using .hPFMD

control=Dict("dt"=>0.01,
            "thermofreq"=>100,
            "zero velocity"=>true,
            "steps"=>1000,
            "velocity verlet"=>true,
            "LAMMPSTrj"=>Dict("file"=>"monoliquid.lammpstrj","freq"=>1000),
            "BerendsenNVT"=>Dict("tau"=>1.0),
            "Canonical field"=>Dict(
            "χ"=>Dict("1"=>[0.0]),
            "κ"=>0.1,
            "uniform mesh"=>true,
            "update"=>1,
            "Lcell"=>1.0),
            "thermo information"=>Dict( "energy"=>true,
                                       "momentum"=>true,
                                       "write"=>true,
                                        "file"=>"monoliquid.log"))


Waters=Dict()
for i in 1:10
    Waters[string(i)]=Dict("atoms"=>Dict("1"=>Dict("type"=>1,"mass"=>1,"coords"=>rand(Float64, 3)*4)))
end

configs=Dict(
    "box"=>[5.0,5.0,5.0],
    "molecules"=>Waters
    )

Simulate(control,configs)
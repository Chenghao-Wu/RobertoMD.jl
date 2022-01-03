
include("/home/zwq2834/development/hPF-MD.jl/hPF-MD.jl/hPF-MD.jl")

using .hPFMD

control=Dict("dt"=>0.01,
            "thermofreq"=>100,
            "zero velocity"=>true,
            "steps"=>10000,
            "velocity verlet"=>true,
            "LangevinNVT"=>Dict("gamma"=>10.0),
            "Canonical field"=>Dict(
            "χ"=>Dict("1"=>[0.0]),
            "κ"=>0.1,
            "uniform mesh"=>true,
            "update"=>100,
            "Lcell"=>1.0))


Waters=Dict()
for i in 1:10000
    Waters[string(i)]=Dict("atoms"=>Dict("1"=>Dict("type"=>1,"mass"=>1,"coords"=>rand(Float64, 3)*10)))
end

configs=Dict(
    "box"=>[10.0,10.0,10.0],
    "molecules"=>Waters
    )

@time Simulate(control,configs)
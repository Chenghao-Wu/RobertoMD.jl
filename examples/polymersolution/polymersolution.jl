
include("/home/zwq2834/development/hPF-MD.jl/hPF-MD.jl/hPF-MD.jl")

using .hPFMD
using Rotations
using JSON




control=Dict("dt"=>0.01,
            "zero velocity"=>true,
            "steps"=>10000,
            "velocity verlet"=>true,
            "LAMMPSTrj"=>Dict("file"=>"polymersolution.lammpstrj","freq"=>1000),
            "restart"=>Dict("JSONRestart"=>true,"file"=>"polymersolution.json","freq"=>1000),
            "BerendsenNVT"=>Dict("tau"=>1.0),
            "harmonic bond"=>Dict(
            "k"=>Dict("1"=>1000.0),
            "l0"=>Dict("1"=>1.0)
            ),
            "Canonical field"=>Dict(
            "χ"=>Dict(  "1"=>[0.0,10.0],
                        "2"=>[10.0,0.0] ),
            "κ"=>0.2,
            "uniform mesh"=>true,
            "update"=>1,
            "Lcell"=>1.0),
            "thermo information"=>Dict( "freq"=>1,
                                        "energy"=>true,
                                       "momentum"=>true,
                                       "write"=>true,
                                        "file"=>"polymersolution.log"),
            "density"=>0.85,
            "bondlength"=>1.0,
            "chain_length"=>20,
            "num_polymers"=>100,
            "num_water"=>20000
            )


boxsize=((control["num_polymers"]*control["chain_length"]+control["num_water"])/control["density"])^(1.0/3)
bondlength=control["bondlength"]
chain_length=control["chain_length"]

polymersolution=Dict()
for polymer_i in 1:control["num_polymers"]
    polymersolution[string(polymer_i)]=Dict()
    polymersolution[string(polymer_i)]["atoms"]=Dict()
    polymersolution[string(polymer_i)]["bonds"]=Dict()
    pos_init=rand(3)*boxsize
    for monomer_i in 1:chain_length
        r = rand(RotMatrix{3})
        q = QuatRotation(r)
        pos_init=pos_init+q*[bondlength,0.0,0.0]
        type_=1
        polymersolution[string(polymer_i)]["atoms"][string(monomer_i)]=Dict("type"=>type_,"mass"=>1,"coords"=>pos_init)
    end
    for bond_i in 1:chain_length-1
        polymersolution[string(polymer_i)]["bonds"][string(bond_i)]=[1,bond_i,bond_i+1]
    end
end

for wateri in control["num_polymers"]+1:control["num_water"]+control["num_polymers"]
    polymersolution[string(wateri)]=Dict("atoms"=>Dict("1"=>Dict("type"=>2,"mass"=>1,"coords"=>rand(Float64, 3)*boxsize)))
end

configs=Dict(
    "box"=>[boxsize,boxsize,boxsize],
    "molecules"=>polymersolution
    )


#control_file=open("control.json","w")
#JSON.print(control_file,control,4)

#=

configuration_file=open("configuration.json","w")
JSON.print(configuration_file,configs,4)

=#

#configs=JSON.parse(open("configuration.json","r"))

Simulate(control,configs)
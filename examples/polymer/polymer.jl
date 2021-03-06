

using RobertoMD
using Rotations
using JSON




control=Dict("dt"=>0.01,
            "zero velocity"=>true,
            "steps"=>10000,
            "velocity verlet"=>true,
            "LAMMPSTrj"=>Dict("file"=>"polymer.lammpstrj","freq"=>1000),
            "restart"=>Dict("JSONRestart"=>true,"file"=>"polymer.json","freq"=>1000),
            "BerendsenNVT"=>Dict("tau"=>1.0),
            "harmonic bond"=>Dict(
            "k"=>Dict("1"=>1000.0),
            "l0"=>Dict("1"=>1.0)
            ),
            "Canonical field"=>Dict(
            "χ"=>Dict("1"=>[0.0]),
            "κ"=>0.2,
            "uniform mesh"=>true,
            "update"=>1,
            "Lcell"=>1.0),
            "thermo information"=>Dict( "freq"=>1,
                                        "energy"=>true,
                                       "momentum"=>true,
                                       "write"=>true,
                                        "file"=>"polymer.log"),
            "density"=>0.85,
            "bondlength"=>1.0,
            "chain_length"=>10,
            "num_polymers"=>500
            )


boxsize=(control["num_polymers"]*control["chain_length"]/control["density"])^(1.0/3)
bondlength=control["bondlength"]
chain_length=control["chain_length"]

polymer=Dict()
for polymer_i in 1:control["num_polymers"]
    polymer[string(polymer_i)]=Dict()
    polymer[string(polymer_i)]["atoms"]=Dict()
    polymer[string(polymer_i)]["bonds"]=Dict()
    pos_init=rand(3)*boxsize
    for monomer_i in 1:chain_length
        r = rand(RotMatrix{3})
        q = QuatRotation(r)
        pos_init=pos_init+q*[bondlength,0.0,0.0]
        polymer[string(polymer_i)]["atoms"][string(monomer_i)]=Dict("type"=>1,"mass"=>1,"coords"=>pos_init)
    end
    for bond_i in 1:chain_length-1
        polymer[string(polymer_i)]["bonds"][string(bond_i)]=[1,bond_i,bond_i+1]
    end

end

configs=Dict(
    "box"=>[boxsize,boxsize,boxsize],
    "molecules"=>polymer
    )


#control_file=open("control.json","w")
#JSON.print(control_file,control,4)

#=

configuration_file=open("configuration.json","w")
JSON.print(configuration_file,configs,4)

=#

#configs=JSON.parse(open("configuration.json","r"))

Simulate(control,configs)
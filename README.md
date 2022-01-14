# RobertoMD

RobertoMD.jl ([JuliaCon2021](https://live.juliacon.org/talk/ECKGDE)): Rome-Berlin-Tokyo Hybrid Particle-field Molecular-Dynamics Simulation

A massively parallel hybrid particle-field molecular dynamics simulation package written in Julia. **It is aimed to become a productive hPF-MD simulator. However, it has not yet been at that stage**. Benchmarks and tests are welcome.

The implemeted functions are limited for now, including:

* Original hPF-MD algorithm of [Milano *et al.*](https://doi.org/10.1063/1.3142103).
    * Particle-in-cell transformation
    * Periodic cubic fields (simulation box must be same in x,y,z directions)
* Bond interactions: Harmonic
* Thermostats: Langevin, Berendsen
* Velocity Verlet integration
* Periodic boundary conditions in a cubic box.
* Particle-decomposition parallelization (MPI.jl)

## Installation

Julia is required, with Julia v1.5 required to get the latest version of RobertoMD. Install Roberto from the Julia REPL. Enter the package mode by pressing ] and run "add RobertoMD".

## Parallelization
The parallelization is implemented via the "particle decomposition" algorithm, taking advantage of the unqiue feature of the field interaction. Details can be found in the [previous article](https://doi.org/10.1002/jcc.22883) [Ref.2].

In short, the field is a function of the local particle density. This collective variable is a slow variable compared to the displacement of the particles. Therefore, it is possible to update the field every $\Delta t$ steps without losing accuracy. In this way, data exchange between different cores is less frequent compared with domain decomposition.

## Usage

Several example systems, e.g., **monoliquid**, **polymer melt**, **copolymer melt**, **polymer in solution**, can be found in Example folder:

This is the hPF-MD simulation of a simple block copolymer system:

```python
using RobertoMD
using Rotations
using JSON

#control parameters
control=Dict(
    "dt"=>0.01,
    "steps"=>1000,
    "velocity verlet"=>true,
    "LAMMPSTrj"=>Dict("file"=>"copolymer.lammpstrj","freq"=>1000),
    "restart"=>Dict("JSONRestart"=>true,"file"=>"copolymer.json","freq"=>1000),
    "BerendsenNVT"=>Dict("tau"=>1.0),
    "harmonic bond"=>Dict(
    "k"=>Dict("1"=>1000.0),
    "l0"=>Dict("1"=>1.0)
    ),
    "Canonical field"=>Dict(
    "χ"=>Dict(  "1"=>[0.0,5.0],
                "2"=>[5.0,0.0] ),
    "κ"=>0.2,
    "uniform mesh"=>true,
    "update"=>1,
    "Lcell"=>1.0),
    "thermo information"=>Dict( "freq"=>1,
                                "energy"=>true,
                                "momentum"=>true,
                                "write"=>true,
                                "file"=>"copolymer.log"),
    "density"=>0.85,
    "bondlength"=>1.0,
    "chain_length"=>20,
    "num_polymers"=>500,
    "zero velocity"=>true,
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
        type_=1
        if monomer_i>div(control["chain_length"],2)
            type_=2
        end
        polymer[string(polymer_i)]["atoms"][string(monomer_i)]=Dict("type"=>type_,"mass"=>1,"coords"=>pos_init)
    end
    for bond_i in 1:chain_length-1
        polymer[string(polymer_i)]["bonds"][string(bond_i)]=[1,bond_i,bond_i+1]
    end

end

configs=Dict(
    "box"=>[boxsize,boxsize,boxsize],
    "molecules"=>polymer
    )



Simulate(control,configs)
```

## References:

1. Milano, G.; Kawakatsu, T. Hybrid Particle-Field Molecular Dynamics Simulations for Dense Polymer Systems. [J. Chem. Phys. 2009, 130 (21), 214106.](https://doi.org/10.1063/1.3142103)
2. Zhao, Y.; Nicola, A. D.; Kawakatsu, T.; Milano, G. Hybrid Particle-Field Molecular Dynamics Simulations: Parallelization and Benchmarks. [J. Comput. Chem. 2012, 33 (8), 868–880.](https://doi.org/10.1002/jcc.22883)
3. Bore, S. L.; Cascella, M. Hamiltonian and Alias-Free Hybrid Particle–Field Molecular Dynamics. [J. Chem. Phys. 2020, 153 (9), 094106.](https://doi.org/10.1063/5.0020733)
4. Wu, Z.; Kalogirou, A.; De Nicola, A.; Milano, G.; Müller‐Plathe, F. Atomistic Hybrid Particle‐field Molecular Dynamics Combined with: Restoring Entangled Dynamics to Simulations of Polymer Melts. [J. Comput. Chem. 2021, 42 (1), 6–18.](https://doi.org/10.1002/jcc.26428)
5. Wu, Z.; Milano, G.; Müller-Plathe, F. Combination of Hybrid Particle-Field Molecular Dynamics and Slip-Springs for the Efficient Simulation of Coarse-Grained Polymer Models: Static and Dynamic Properties of Polystyrene Melts. [J. Chem. Theory Comput. 2021 17 (1), 474-487](https://doi.org/10.1021/acs.jctc.0c00954)
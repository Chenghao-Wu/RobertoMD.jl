# RobertoMD

RobertoMD.jl: Rome-Berlin-Tokyo Hybrid Particle-field Molecular-Dynamics Simulation

A massively parallel hybrid particle-field molecular dynamics simulation package written in Julia

The parallelization is implemented via the "particle decomposition" algorithm, taking advantage of the unqiue feature of the field interaction. Details can be found in the [previous article](https://doi.org/10.1002/jcc.22883) [Ref.2].

In short, the field is a function of the local particle density. This collective variable is a slow variable compared to the displacement of the particles. Therefore, it is possible to update the field every $\Delta t$ steps without losing accuracy. In this way, data exchange between different cores is less frequent compared with domain decomposition.

References:

1. [Milano, G.; Kawakatsu, T. Hybrid Particle-Field Molecular Dynamics Simulations for Dense Polymer Systems. J. Chem. Phys. 2009, 130 (21), 214106.](https://doi.org/10.1063/1.3142103)
2. [Zhao, Y.; Nicola, A. D.; Kawakatsu, T.; Milano, G. Hybrid Particle-Field Molecular Dynamics Simulations: Parallelization and Benchmarks. J. Comput. Chem. 2012, 33 (8), 868–880.](https://doi.org/10.1002/jcc.22883)
3. [Bore, S. L.; Cascella, M. Hamiltonian and Alias-Free Hybrid Particle–Field Molecular Dynamics. J. Chem. Phys. 2020, 153 (9), 094106.](https://doi.org/10.1063/5.0020733)
4. [Wu, Z.; Kalogirou, A.; De Nicola, A.; Milano, G.; Müller‐Plathe, F. Atomistic Hybrid Particle‐field Molecular Dynamics Combined with: Restoring Entangled Dynamics to Simulations of Polymer Melts. J. Comput. Chem. 2021, 42 (1), 6–18.](https://doi.org/10.1002/jcc.26428)
5. [Wu, Z.; Milano, G.; Müller-Plathe, F. Combination of Hybrid Particle-Field Molecular Dynamics and Slip-Springs for the Efficient Simulation of Coarse-Grained Polymer Models: Static and Dynamic Properties of Polystyrene Melts. J. Chem. Theory Comput. 2020.](https://doi.org/10.1021/acs.jctc.0c00954)
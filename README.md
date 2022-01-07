# RobertoMD

RobertoMD.jl: Rome-Berlin-Tokyo Hybrid Particle-field Molecular-Dynamics Simulation

A massively parallel hybrid particle-field molecular dynamics simulation package written in Julia

The parallelization is implemented via the "particle decomposition" algorithm, taking advantage of the unqiue feature of the field interaction. Details can be found in the previous article.

In short, the field is a function of the local particle density. This collective variable is a slow variable compared to the displacement of the particles. Therefore, it is possible to update the field every $\Delta t$ steps without losing accuracy. In this way, data exchange between different cores is less frequent compared with domain decomposition.

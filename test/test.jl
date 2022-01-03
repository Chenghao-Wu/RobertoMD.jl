
include("/home/zwq2834/development/hPF-MD.jl/hPF-MD.jl/hPF-MD.jl")

using .hPFMD

inputs=ReadInput("input.json")
configs=ReadConfig("configuration.json")

@time Simulate(inputs,configs)
module SQRLab

using IonSim
using Distributions
using LsqFit
using Plots

include("Trap.jl")
include("Simulation.jl")

using .Trap
using .Simulation

export create_standard_chamber, run_simulation

end
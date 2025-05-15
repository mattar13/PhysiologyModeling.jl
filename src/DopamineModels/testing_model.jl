using DifferentialEquations
using Pkg; Pkg.activate("test")
using GLMakie

include("parameters.jl")
include("auxillary_functions.jl")
include("models.jl")

# ---------- simulation ----------
u0 = [-70.0, 0.05, 0.0, 0.0, 0.0, 2.0, 0.2]   # initial conditions
tspan = (0.0, 500.0)                          # ms
prob  = ODEProblem(dopaminergic_autoreceptor!, u0, tspan, p)
sol   = solve(prob, Tsit5(); saveat=0.5)      # 0.5 ms sampling
using Revise
using PhysiologyModeling
using PhysiologyPlotting
using PyPlot
#Load the dimensions at the top of the stack
@parameters t x y

include("..\\src\\StarburstModel\\Variables.jl")
include("..\\src\\StarburstModel\\Equations.jl")
tspan = (0.0, 100.0)
@named SAC_ODE = ODESystem(SAC_eqs, t, states, parameters)
@named SAC_SDE = SDESystem(SAC_ODE, SAC_noise_eqs, tspan = (0.0, 60e3))
prob = SDEProblem(SAC_SDE)
sol = solve(prob, SOSRI())

#Make an updated Experiment
data = Experiment(sol)
plot_experiment(data, sweeps = [2, 11])
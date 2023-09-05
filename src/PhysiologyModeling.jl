module PhysiologyModeling

# Write your package code here.
using DifferentialEquations, ModelingToolkit
#These exports are used in modelling
export @variables, @parameters, Differential
export @connector, @component, @unpack
export Flow, Equation, compose, extend, connect
export ODAEProblem, ODEProblem, Tsit5
export @named, ODESystem
export structural_simplify, ODEProblem, solve, remake 
export PresetTimeCallback 
#Export some commonly used algorithims
export Rodas5

@variables t #Define the time variables which will stream everything in the model

include("InitialConditions/conditions.jl")
include("InitialConditions/parameters.jl")
include("Components/auxillary_functions.jl") #This loads all of the components
include("Components/components.jl") #This loads all of the components
export t, LeakySome, Cable

#This section deals with parameters and contions
#Eventually PhysiologyPlotting will include some things we need to plot everything
end

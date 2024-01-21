module PhysiologyModeling

using Requires
using ElectroPhysiology
#using PhysiologyPlotting
# Write your package code here.
using LinearAlgebra,SparseArrays
using DifferentialEquations

using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

#Export some commonly used algorithims
export Tsit5, Rodas5, ROS3P, TRBDF2, KenCarp47
export AutoTsit5, Rosenbrock23
export SOSRI, SOSRA, SOSRA2
export solve
export SDEProblem, ODEProblem

include("utilities.jl")
export Experiment

include("StarburstModel/AuxillaryFunctions.jl")

include("StarburstModel/Variables.jl")
export SAC_u0_dict, SAC_p0_dict
export keys_u0, keys_p0
export vals_u0, vals_p0
export nt_u0, nt_p0
export extract_p0, extract_u0, extract_dict

include("StarburstModel/Models.jl")
export SAC_ODE, SAC_ODE_NT_CLAMP, SAC_ODE_STIM
export ∇α, DIFFUSION_MODEL, DIFFUSION_NOISE
export SAC_PDE, SAC_PDE_STIM 
export SAC_GAP #Make a model for a gap junction
export SAC_ODE_Compartment
export noise1D, noise2D

include("StarburstModel/Mapping.jl")
export even_map
export CellMap
export ring_circle_overlap_area
export make_GPU

#This section deals with parameters and contions
#Eventually PhysiologyPlotting will include some things we need to plot everything

function __init__()
    @require CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba" begin
        println("GPU mode enabled")
        using .CUDA
        include("make_GPU.jl")
        export make_GPU, CellMap_GPU, SAC_PDE_GPU
    end
end


end

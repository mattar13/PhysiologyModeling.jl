module PhysiologyModeling

using Requires
using ElectroPhysiology
#using PhysiologyPlotting
# Write your package code here.
using LinearAlgebra, SparseArrays
using NLsolve
#using DiffEqCallbacks #This is necessary for inserting conductances into the model

using Logging: global_logger
using TerminalLoggers: TerminalLogger

global_logger(TerminalLogger())

#Export some commonly used algorithims
export EM, Euler
export Tsit5, Rodas5, ROS3P, TRBDF2, KenCarp47
export AutoTsit5, Rosenbrock23
export remake
export EnsembleProblem, EnsembleDistributed
export SOSRI, SOSRA, SOSRA2
export solve
export SDEProblem, ODEProblem

include("utilities.jl")
export Experiment
export extract_p0, extract_u0, extract_dict

include("StarburstModel/AuxillaryFunctions.jl")
export gauss_pulse #generating a pulse of glutamate

include("StarburstModel/Variables.jl")
export SAC_u0_dict, SAC_p0_dict
export extract_u0, extract_p0
export par_idx

include("StarburstModel/Models.jl")
export SAC_ODE, SAC_ODE_NT_CLAMP
export SAC_ODE_IC, SAC_ODE_VC
export ∇α, DIFFUSION_MODEL, DIFFUSION_NOISE
export SAC_PDE, SAC_PDE_STIM 
export SAC_GAP #Make a model for a gap junction
export SAC_ODE_Compartment
export noise1D, noise2D
export DynamicSAC
export apply_glutamate_affect!

include("StarburstModel/Mapping.jl")
export CellMap

export create_even_map, create_random_map
export create_ring_map, create_dendrogram_map

export connect_k_neighbors, connect_neighbors_radius
export connection_matrix

export RING, CIRC, RING_CIRC, LIN
export RING_CIRC_BIAS

export find_angle, calculate_exponential_bias

export make_GPU

#Eventually PhysiologyPlotting will include some things we need to plot everything

function __init__()
    @require CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba" begin
        println("GPU mode enabled")
        using .CUDA
        include("make_GPU.jl")
        export make_GPU, CellMap_GPU, SAC_PDE_GPU 
        export DIFFUSION_MODEL_GPU
    end

    #put the exports under here
    @require BifurcationKit = "0f109fa4-8a5d-4b75-95aa-f515264e7665" begin
        println("BifurcationKit loaded")
        
        include("DynamicalAnalysis/ensemble_functions.jl")
        export OneVarEnsembleProb, InitialCond_Param_EnsembleProb
        
        include("DynamicalAnalysis/phase_plane_analysis.jl")
        export phase_plane, find_nullclines
        
        include("DynamicalAnalysis/equilibria_analysis.jl")
        export find_fixed_points, find_equilibria
        @require DynamicalSystems = "61744808-ddfa-5f27-97ff-6e42cc95d634" begin
            println("Dynamical functions loaded")
        end
    end
end


end

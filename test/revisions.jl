using Revise
using PhysiologyModeling
using Pkg; Pkg.activate("test") #Activate the testing environment
#I think a good next step is to add in the ability to open and compare data versus physiologyical data
# create a model output here
using ElectroPhysiology
using PhysiologyPlotting, GLMakie

f = Figure(size = (2.5, 2.5))
ax1 = Axis(f[1,1], 
     title = "Experiment Plot Test",
     xlabel = "Time (ms)", 
     ylabel = "Response (μV)"
)
#First try opening some data
#find a good data point
folder = raw"F:\SAC_project_backup_Feller_Lab\SAC-project_ephys_data\01132023-P10\B2-KO-nGFP\cell1_SAC"
file = "gap_free_voltage_0000.abf"
path = joinpath(folder, file)
data = readABF(path, channels = ["I_MTest 1"], stimulus_name = nothing)
size(data)
data.t

downsample!(data, 10.0)
plot_experiment(ax1, data)

display(f)
#1) Run the model using the default settingstspan = (0.0, 10e3)
p0 = extract_p0(SAC_p0_dict)
prob = SDEProblem(SAC_ODE, noise1D, vals_u0, (0.0, 300e3), p0)
@time sol = solve(prob, SOSRI(), reltol = 0.01, abstol = 0.01, progress = true, progress_steps = 1)
sim_exp = Experiment(sol) #Turn this into an experiment

f = Figure(size = (2.5, 2.5))
ax1 = Axis(f[1,1], 
     title = "Experiment Plot Test",
     xlabel = "Time (ms)", 
     ylabel = "Response (μV)"
)
plot_experiment(ax1, sim_exp)


#%% This section is testing the GPU 
#1) Run the model using the default settings
#I want to use GPUs to improve the speed of the diffusion reaction. 

using Revise, Profile, ProfileSVG
using CUDA 
using BenchmarkTools
using PhysiologyModeling
using PhysiologyPlotting
using GLMakie
using SparseArrays, LinearAlgebra
import PhysiologyModeling: Φe

#1) determine the domains and spacing of cells. 
domain_x = (xmin, xmax) = (0.0, 10.0)
domain_y = (ymin, ymax) = (0.0, 10.0)
dx = dy = 0.05 #Mean distribution is 40-50 micron (WR taylor et al)

#2) create the map of cells and their radii
cells = even_map(xmin = xmin, dx = dx, xmax = xmax, ymin = ymin, dy = dy, ymax = ymax)
radii = fill(0.200, size(cells, 1)) #Switch this on to get constant radii
cell_map = CellMap(cells, radii);

#3) Define the initial state and timespan
u = ones(size(cell_map.connections, 1))
du = similar(u)

uG = ones(size(cell_map.connections, 1)) |> CuArray
duG = similar(uG)
cell_connect_G = cell_map.connections |> CuSparseMatrixCSC

function testf(du, u, S)
    mul!(du, S, u);
end

function testf_GPU(du, u, S)
    #CUDA.@sync mul!(du, S, u);
    CUDA.@sync du .= S * u
end

@time cell_map.connections * du;

@time testf(du, u, cell_map.connections);
@time testf_GPU(duG, uG, cell_connect_G);
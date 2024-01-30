using Revise
using ElectroPhysiology
using PhysiologyModeling
using Pkg; Pkg.activate("test") #Activate the testing environment
using PhysiologyPlotting, GLMakie
using BenchmarkTools

#%%
tspan = (0.0, 4000.0)
SAC_p0_dict["I_app"] = 5.0
#%% Section A: Testing models against physiology data
#I think a good next step is to add in the ability to open and compare data versus physiologyical data
root = raw"E:\Data\Patching"
file = "2024_01_25_ChAT-RFP_DSGC/Cell1/24125000.abf"
filename = joinpath(root, file)
data = readABF(filename)
import ElectroPhysiology.create_signal_waveform!
create_signal_waveform!(data, "Analog 0")


f = Figure(size = (200, 200))
ax1 = Axis(f[1,1], 
     title = "Experiment Plot Test",
     xlabel = "Time (ms)", 
     ylabel = "Response (mV)"
)
ax2 = Axis(f[1,2])

plot_experiment(ax1, exp; channels = 1)
plot_experiment(ax1, exp; channels = 2)

#First try opening some data
#find a good data point
data = readABF(filename)
size(data)


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

# Step 1. Set up all parameters for the ODE
tspan = (0.0, 120e3)

SAC_p0_dict["I_app"] = 0.0
SAC_p0_dict["g_GABA"] = 0.0
SAC_p0_dict["g_ACh"] = 0.0
SAC_p0_dict["g_W"] = 0.0
SAC_p0_dict["g_TREK"] = 0.0
p0 = extract_p0(SAC_p0_dict)
u0 = extract_u0(SAC_u0_dict)
prob = SDEProblem(SAC_ODE_STIM, noise1D, u0, tspan, p0)

xlims = (-90.0, 10.0)
ylims = (-0.10, 5.0)
xmap = LinRange(xlims[1], xlims[2], 100)
ymap = LinRange(ylims[1], ylims[2], 100)
@time zplane = phase_plane(prob, xmap, ymap);
#arrows(xmap, ymap, zplane[:,:,1], zplane[:,:,2], align = :center, arrowsize = 1, lengthscale = 0.03,)

prob |> typeof |> fieldnames

find_equilibria(prob, xlims, ylims, verbose = true)

#%% Section B: Running using GPUs and phys data
using CUDA
using SparseArrays
using LinearAlgebra

#1) determine the domains and spacing of cells. 
domain_x = (xmin, xmax) = (0.0, 10.0)
domain_y = (ymin, ymax) = (0.0, 10.0)
dx = dy = 0.05 #Mean distribution is 40-50 micron (WR taylor et al)

#2) create the map of cells and their radii
cells = even_map(xmin = xmin, dx = dx, xmax = xmax, ymin = ymin, dy = dy, ymax = ymax)
radii = fill(0.200, size(cells, 1)) #Switch this on to get constant radii
cell_map = make_GPU(CellMap(cells, radii));
u0 = vcat(fill(vals_u0, size(cells, 1))'...) |> CuArray{Float32} #Generate a new initial conditions

tspan = (0.0, 10000.0)
f_PDE(du, u, p, t) = SAC_PDE_GPU(du, u, p, t, cell_map)
prob = SDEProblem(f_PDE, noise2D, u0, tspan, p0)
@Time sol = solve(prob, SOSRI(), reltol = 2e-2, abstol = 2e-2, progress=true, progress_steps=1);

#%% Section A: Testing models against physiology data
#I think a good next step is to add in the ability to open and compare data versus physiologyical data

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
downsample!(data, 10.0)
plot_experiment(ax1, data)

display(f)
#1) Run the model using the default settingstspan = (0.0, 10e3)
p0 = extract_p0(SAC_p0_dict)
prob = SDEProblem(SAC_ODE, noise1D, vals_u0, (0.0, 300e3), p0)
@Time sol = solve(prob, SOSRI(), reltol = 0.01, abstol = 0.01, progress = true, progress_steps = 1)
sim_exp = Experiment(sol) #Turn this into an experiment

f = Figure(size = (2.5, 2.5))
ax1 = Axis(f[1,1], 
     title = "Experiment Plot Test",
     xlabel = "Time (ms)", 
     ylabel = "Response (μV)"
)
plot_experiment(ax1, sim_exp)

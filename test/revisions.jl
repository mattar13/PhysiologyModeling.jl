using Revise
using ElectroPhysiology
using PhysiologyModeling
PhysiologyModeling.__init__()
using Pkg; Pkg.activate("test") #Activate the testing environment
using PhysiologyPlotting, GLMakie
using BenchmarkTools

#3 Goals
#GOAL) Analysis of Physiological data
#GOAL) Speeding up the model
#GOAL) Analyzing the model
#GOAL) Testing models against physiology data

#%% GOAL: Testing models against physiology data
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


#%% GOAL: Speeding up the model
using CUDA
CUDA.allowscalar(false)
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
SAC_p0_dict["g_GABA"] = 0.0
p0 = extract_p0(SAC_p0_dict)

#3) Define the problem
tspan = (0.0, 100.0)
f_PDE(du, u, p, t) = SAC_PDE_GPU(du, u, p, t, cell_map)
prob = SDEProblem(f_PDE, noise2D, u0, tspan, p0)

#%% Pause here before running the model
@time sol = solve(prob, SOSRI(), reltol = 2e-3, abstol= 2e-3, progress=true, progress_steps=1)
sol.t

#%% Plot the solutions
CUDA.allowscalar(true)
fSDE = Figure(size = (1800, 800))
ax1 = Axis(fSDE[1,1], title = "Voltage (Vt)")
ax2 = Axis(fSDE[2,1], title = "K Repol. (Nt)")
ax3 = Axis(fSDE[3,1], title = "Na Gating (Mt)")
ax4 = Axis(fSDE[4,1], title = "Na Close (Ht)")

ax5 = Axis(fSDE[1,2], title = "Calcium (Ct)")
ax6 = Axis(fSDE[2,2], title = "cAMP (At)")
ax7 = Axis(fSDE[3,2], title = "TREK (Bt)")

ax8 = Axis(fSDE[1,3], title = "ACh (Et)")
ax9 = Axis(fSDE[2,3], title = "GABA (It)")

ax10 = Axis(fSDE[1,4], title = "Noise (Wt)")
Time = LinRange(sol.t[1], sol.t[end], 1000)
for i in rand(1:size(sol, 1), 100)
     lines!(ax1, Time, map(t -> sol(t)[i, 1], Time))
     lines!(ax2, Time, map(t -> sol(t)[i, 2], Time))
     lines!(ax3, Time, map(t -> sol(t)[i, 3], Time))
     lines!(ax4, Time, map(t -> sol(t)[i, 4], Time))
     lines!(ax5, Time, map(t -> sol(t)[i, 5], Time))
     lines!(ax6, Time, map(t -> sol(t)[i, 6], Time))
     lines!(ax7, Time, map(t -> sol(t)[i, 7], Time))
     lines!(ax8, Time, map(t -> sol(t)[i, 8], Time))
     lines!(ax9, Time, map(t -> sol(t)[i, 9], Time))
     lines!(ax10, Time, map(t -> sol(t)[i, 10], Time))
end
display(fSDE)

#%% GOAL:  Analyzing the model
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
find_equilibria(prob, xlims, ylims, verbose = true)
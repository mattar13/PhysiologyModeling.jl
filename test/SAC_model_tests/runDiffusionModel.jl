using Revise, BenchmarkTools
using PhysiologyModeling
using PhysiologyPlotting
using GLMakie
using SparseArrays, LinearAlgebra
import PhysiologyModeling: Φe, find_neighbors_radius, create_sparse_matrix
import PhysiologyModeling:OLD_∇α
#%% Set up the PDE
#Here is all of the ways to set up morphology
domain_x = (xmin, xmax) = (0.0, 1.0)
domain_y = (ymin, ymax) = (0.0, 1.0)
dx = dy = 0.05 #Mean distribution is 40-50 micron (WR taylor et al)

cells = even_map(xmin = xmin, dx = dx, xmax = xmax, ymin = ymin, dy = dy, ymax = ymax)
#radii = rand(0.10:0.01:0.20, size(cells, 1))
radii = fill(0.200, size(cells, 1))
cell_map = CellMap(cells, radii, max_strength = 0.005);

#Define the initial state and timespan
u0 = zeros(size(cell_map.connections, 1))
mid = round(Int64, size(cell_map.xs, 1)/2)+1

cell_map.xs

#Run model
tspan = (0.0, 5000)
probSDE = SDEProblem(DIFFUSION_MODEL, DIFFUSION_NOISE, u0, tspan, cell_map)
@time sol = solve(probSDE, SOSRI(), reltol=0.01, abstol=0.01, progress=true, progress_steps=1)

#%% Animate the soluiton
fDIFF = Figure(resolution = (1000,1000))
ax1 = Axis3(fDIFF[1,1]; aspect=(1, 1, 1), title = "T = 0.00")

surf2 = surface!(ax1, cells[:, 1], cells[:, 2], sol(0.0), colorrange = (0.0, 0.5), markersize = 35.0)

xlims!(ax1, (xmin, xmax))
ylims!(ax1, (ymin, ymax))
zlims!(ax1, (0.0, maximum(sol)))
n_frames = 1000
animate_t = LinRange(0.0, sol.t[end], n_frames)
dt = animate_t[2] - animate_t[1]
fps = round(Int64, (1/dt) * 1000)
record(fDIFF, "test/SAC_model_tests/diffusion_animation.mp4", animate_t, framerate = fps) do t
	println(t)
	u = sol(t)
	ax1.title = "T = $(round(t/1000, digits = 3))s"
	surf2[3] = u
end
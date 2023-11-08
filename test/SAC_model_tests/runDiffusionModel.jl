using Revise
using PhysiologyModeling
using PhysiologyPlotting
using GLMakie
using SparseArrays
import PhysiologyModeling.Φe
#%% Set up the PDE
#Here is all of the ways to set up morphology
dpi = 400
xmin = ymin = 0.0
xmax = ymax = 0.7
dx = dy = 0.05 #Mean distribution is 40-50 micron (WR taylor et al)
cells = even_map(xmin = xmin, dx = dx, xmax = xmax, ymin = ymin, dy = dy, ymax = ymax)
#radii = rand(0.1:0.01:0.20, size(cells, 1))
radii = fill(0.200, size(cells, 1))
cell_map = CellMap(cells, radii, max_strength = 0.005);

#Define the initial state and timespan
u0 = zeros(size(cell_map.connections, 1))
mid = round(Int64, size(cell_map.xs, 1)/2)+1
u0[mid] = 1.0 #Add a "bolus" of material at the center

tspan = (0.0, 100000.0)
#prob = ODEProblem(DIFFUSION_MODEL, u0, tspan, cell_map)
#make a SDE problem
probSDE = SDEProblem(DIFFUSION_MODEL, DIFFUSION_NOISE, u0, tspan, cell_map)
# Solve the differential equation
#@time sol = solve(prob, TRBDF2(), reltol=1e-3, abstol=0.01, progress=true, progress_steps=1)
@time sol = solve(probSDE, SOSRI(), reltol=1e-3, abstol=0.01, progress=true, progress_steps=1)

sol |> maximum

#%% Plot the figure
# Animation settings
n_frames = 1000
animate_t = LinRange(0.0, sol.t[end], n_frames)
dt = animate_t[2] - animate_t[1]
fps = round(Int64, (1/dt) * 1000)

fDIFF = Figure(resolution = (1000,1000))
ax11 = Axis3(fDIFF[1,1]; aspect=(1, 1, 1))
ax12 = Axis3(fDIFF[1,2]; aspect=(1, 1, 1))
#ax13 = Axis3(fDIFF[1,3]; aspect=(1, 1, 1))

ax21 = Axis(fDIFF[2,1])
ax22 = Axis(fDIFF[2,2])
#ax23 = Axis(fDIFF[2,3])

surface!(ax11, cell_map.xs, cell_map.ys, sol(0.0), colorrange = (0.0, 0.5))
scatter!(ax21, cell_map.xs, cell_map.ys, color = sol(0.0), colorrange = (0.0, 0.5),markersize = 15.0)
scatter!(ax21, cell_map.xs[mid], cell_map.ys[mid], strokewidth = 1, color = :transparent, markersize = cell_map.radius[mid]*2, markerspace = :data)
surf2 = surface!(ax12, cells[:, 1], cells[:, 2], sol(0.0), colorrange = (0.0, 0.5), markersize = 35.0)
scat2 = scatter!(ax22, cells[:, 1], cells[:, 2], color = sol(0.0), colorrange = (0.0, 0.5), markersize = 15.0)
scatter!(ax22, cell_map.xs[mid], cell_map.ys[mid], strokewidth = 1, color = :transparent, markersize = cell_map.radius[mid]*2, markerspace = :data)
xlims!(ax11, (xmin, xmax))
ylims!(ax11, (ymin, ymax))
zlims!(ax11, (0.0, 1.0))

xlims!(ax12, (xmin, xmax))
ylims!(ax12, (ymin, ymax))
zlims!(ax12, (0.0, 1.0))

record(fDIFF, "test/SAC_model_tests/diffusion_animation.mp4", animate_t, framerate = fps) do t
	#println(t)
	u = sol(t)
	surf2[3] = u
	scat2.color = u
end
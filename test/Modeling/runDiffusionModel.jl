using Revise
using PhysiologyModeling
using Pkg; Pkg.activate("test")

using CUDA
using PhysiologyPlotting
using GLMakie
using SparseArrays, LinearAlgebra

#%%
#1) determine the domains and spacing of cells. 
domain_x = (xmin, xmax) = (-1.0, 1.0)
domain_y = (ymin, ymax) = (-1.0, 1.0)
dx = dy = 0.05 #Mean distribution is 40-50 micron (WR taylor et al)

#2) create the map of cells and their radii
xs, ys = create_even_map(xmin = xmin, dx = dx, xmax = xmax, ymin = ymin, dy = dy, ymax = ymax)
length(xs)
connection_list = connect_neighbors_radius(xs, ys, 0.18, self_connecting = false)
connections = connection_matrix(connection_list, m = length(xs), n = length(ys))
dist_func(d) = ring_circle_overlap_area(d; density = 10.0, r_inner = 0.1, r_outer = 0.15, r_circle = 0.18)
cell_map = CellMap(xs, ys, connections; distance_function = dist_func) |> make_GPU;

#3) Define the initial state and timespan
u0 = zeros(size(cell_map.connections, 1)) |> CuArray{Float32}
mid = ceil(Int64, length(xs)/2)

u0[mid] =  1.0
#4) Run model
tspan = (0.0, 1.0)
f_diffuse(du, u, p, t) = DIFFUSION_MODEL_GPU(du, u, p, t; active_cell = mid)
probSDE = SDEProblem(f_diffuse, DIFFUSION_NOISE, u0, tspan, cell_map)
@time sol = solve(probSDE, SOSRI(), reltol=1e-1, abstol=1e-1, progress=true, progress_steps=1)

# Run figure 2, the diffusion animation
fig2 = Figure(size = (1000,1000))
ax1 = Axis3(fig2[1,1]; aspect=(1, 1, 1), title = "T = 0.00")

surf2 = surface!(ax1, xs, ys, sol(0.0), colorrange = (0.0, 0.1), markersize = 35.0)

xlims!(ax1, (xmin, xmax))
ylims!(ax1, (ymin, ymax))
zlims!(ax1, (0.0, 0.1))
n_frames = 1000
animate_t = LinRange(0.0, sol.t[end], n_frames)
dt = animate_t[2] - animate_t[1]
fps = (1/dt)*5
GLMakie.record(fig2, "test/Modeling/Results/diffusion_animation.mp4", animate_t, framerate = fps) do t
	println(t)
	u = sol(t)
	ax1.title = "T = $(round(t, digits = 3))s"
	surf2[3] = u
end

#%% Figure 1 showing the map
fig1 = Figure(size = (500,500))
ax1a = Axis(fig1[1,1], aspect = 1)

for (r,c,d) in connection_list
    lines!(ax1a, [xs[r], xs[c]], [ys[r], ys[c]], color = [d], colorrange = (minimum(data), maximum(data)), alpha = 0.8)
end
texts = map(i -> "$i", length(xs))
n_stops = 4
x_stops = LinRange(minimum(xs), maximum(xs), n_stops)

for i in 1:n_stops-1
    sctV = scatter!(ax1a, xs[findall(x_stops[i] .<= xs .<= x_stops[i+1])], ys[findall(x_stops[i] .<= xs .<= x_stops[i+1])], 
        markersize = 35.0, color = [i], colorrange = (1, n_stops), depth = 2 
    )
end
sctV = scatter!(ax1a, xs, ys, markersize = 25.0, color = :gold, depth = 2)

text!(ax1a, xs, ys, text = texts, align = (:center, :center))
display(fig1)
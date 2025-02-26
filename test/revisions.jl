using Pkg;Pkg.activate(".")
using Revise
using ElectroPhysiology, PhysiologyModeling
using DiffEqCallbacks
using Pkg;Pkg.activate("test")
using GLMakie, PhysiologyPlotting
using CUDA

#%%create a map from preselected dictionary of points
# 1) determine the domains and spacing of cells. 
domain_x = (xmin, xmax) = (0.0, 0.0975) #This is a simulation for a retina 5mm in diameter
domain_y = (ymin, ymax) = (0.0, 0.0975)
dx = dy = 0.02 #Mean distribution is 40-50 micron (WR taylor et al)

#2) create a random distribution of cells and their radii
#The density of SACs in the retina is around 1200 per mm2. So if we have 5mm2 1200 * 5 = 6000
n_cells = 10   #Really pushing the model
xs, ys = create_random_map(n_cells, 
     xmin = xmin, dx = dx, xmax = xmax, 
     ymin = ymin, dy = dy, ymax = ymax
)
connection_list = connect_neighbors_radius(xs, ys, 0.08)
connections = connection_matrix(connection_list, m = length(xs), n = length(ys))

#Make the map for ACh
cell_map_GAP = CellMap(xs, ys, connections; distance_function = (x, y) -> -0.5) |> make_GPU;

#Make the map for ACh
ACH_dist_func(p1, p2) = RING_CIRC(p1, p2) 
cell_map_ACH = CellMap(xs, ys, connections; distance_function = ACH_dist_func) |> make_GPU;
#Make the map for GABA
bias_func(p1, p2) = calculate_exponential_bias(find_angle(p1, p2), 90.0)
GABA_dist_func(p1, p2) = RING_CIRC(p1, p2; density = 0.05) * bias_func(p1, p2)
cell_map_GABA = CellMap(xs, ys, connections; distance_function = GABA_dist_func) |> make_GPU;

#[run the model]____________________________________________________________________________#
#Load parameters
p0_dict = SAC_p0_dict()
p0_dict["I_app"] = 4.0
p0_dict["g_ACh"] = 0.0
p0_dict["g_GABA"] = 0.0
p0_dict["gGAP"] = 10.0
p0 = extract_p0(p0_dict) 

#Load initial conditions
u0_dict= SAC_u0_dict(mode = :PDE, n_cells = n_cells)
u0 = extract_u0(u0_dict) |> CuArray{Float32} 
#3) Define the problem
tspan = (0.0, 100e3)

dosetimes = collect(5e3:1000:7e3)
affect!(integrator) = integrator.u[2, 2] = 10.0
cb = PresetTimeCallback(dosetimes, affect!)

f_PDE(du, u, p, t) = SAC_PDE(du, u, p, t, cell_map_GAP, cell_map_ACH, cell_map_GABA) #for now diffusion is the same in both directions
prob = SDEProblem(f_PDE, noise2D, u0, tspan, p0)
@time sol = solve(prob, 
     SOSRA(),
     reltol =  0.1, abstol = 0.1, 
     #tstops = dosetimes, callback = cb, 
     progress=true, progress_steps=1
)
#save("data.jld", "initial_cond", sol[end])
 
#%%====================================[Plot the solution]====================================#
using Dates
CUDA.allowscalar(true) #allow GPU operations to be offloaded to CPU 
start_time = sol.t[1]
end_time = sol.t[end]
Time = start_time:10:end_time
It = hcat(map(t -> sol(t)[:,1], Time)...)|>Array
vt = hcat(map(t -> sol(t)[:,2], Time)...)|>Array
nt = hcat(map(t -> sol(t)[:,3], Time)...)|>Array
mt = hcat(map(t -> sol(t)[:,4], Time)...)|>Array
ht = hcat(map(t -> sol(t)[:,5], Time)...)|>Array
ct = hcat(map(t -> sol(t)[:,6], Time)...)|>Array
at = hcat(map(t -> sol(t)[:,7], Time)...)|>Array
bt = hcat(map(t -> sol(t)[:,8], Time)...)|>Array
et = hcat(map(t -> sol(t)[:,9], Time)...)|>Array
it = hcat(map(t -> sol(t)[:,10], Time)...)|>Array
gt = hcat(map(t -> sol(t)[:,11], Time)...)|>Array
dt = hcat(map(t -> sol(t)[:,12], Time)...)|>Array
qt = hcat(map(t -> sol(t)[:,13], Time)...)|>Array
Wt = hcat(map(t -> sol(t)[:,14], Time)...)|>Array


fig1 = Figure(size = (1800, 800))
ax1a = Axis(fig1[1,1], title = "I_ext (pA)")
ax2a = Axis(fig1[2,1], title = "Voltage (Vt)")
ax3a = Axis(fig1[3,1], title = "Noise (Wt)")

ax1b = Axis(fig1[1,2], title = "K Repol. (Nt)")
ax2b = Axis(fig1[2,2], title = "Na Gating (Mt)")
ax3b = Axis(fig1[3,2], title = "Na Close (Ht)")

ax1c = Axis(fig1[1,3], title = "Calcium (Ct)")
ax2c = Axis(fig1[2,3], title = "cAMP (At)")
ax3c = Axis(fig1[3,3], title = "TREK (Bt)")

ax1d = Axis(fig1[1,4], title = "ACh (Et)")
ax2d = Axis(fig1[2,4], title = "GABA (It)")
ax3d = Axis(fig1[3,4], title = "Glutamate (Gt)")

for n in 1:n_cells
     lines!(ax1a, Time, It[n, :])
     lines!(ax2a, Time, vt[n, :])
     lines!(ax3a, Time, Wt[n, :])
     
     lines!(ax1b, Time, nt[n, :])
     lines!(ax2b, Time, mt[n, :])
     lines!(ax3b, Time, ht[n, :])
     
     lines!(ax1c, Time, ct[n, :])
     lines!(ax2c, Time, at[n, :])
     lines!(ax3c, Time, bt[n, :])
     
     lines!(ax1d, Time, et[n, :])
     lines!(ax2d, Time, it[n, :])
     lines!(ax3d, Time, gt[n, :])
end
display(fig1)
save_fn1 = "test/$(Dates.today())_MODEL_TRACES.png"
save(save_fn1, fig1)

#%% [Create the cell maps]__________________________________________________________________________________#
fig2 = Figure(size = (400, 800))
ax1 = Axis(fig2[1,1])
ax2 = Axis(fig2[2,1])
sctV = scatter!(ax1, xs, ys, color = sol(start_time)[:,6]|>Array, colorrange = (0.0, maximum(ct)), markersize = 20.0)
for n in 1:n_cells
     lines!(ax2, Time, ct[n, :])
end
ticker = vlines!([0.0])

n_frames = 1000
animate_t = LinRange(start_time, end_time, n_frames)
dt = animate_t[2] - animate_t[1]
fps = round(Int64, (1/dt) * 1000)

save_fn2 = "test/$(Dates.today())_MODEL_VIDEO.mp4"
GLMakie.record(fig2, save_fn2, animate_t, framerate = fps) do t
	println(t)
	c = Array{Float64}(sol(t)[:, 6]) .|> Float64
	sctV.color = c
     ticker[1] = [t]
end
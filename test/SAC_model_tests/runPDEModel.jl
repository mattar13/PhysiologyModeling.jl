using Revise
using ElectroPhysiology
using PhysiologyModeling
import PhysiologyModeling: CVODE_BDF, ring
using Pkg; Pkg.activate("test") #Activate the testing environment
using PhysiologyPlotting, GLMakie
using CUDA
CUDA.allowscalar(false)
using SparseArrays
using LinearAlgebra

#%%

#%% 1) determine the domains and spacing of cells. 
domain_x = (xmin, xmax) = (0.0, 1.0)
domain_y = (ymin, ymax) = (0.0, 1.0)
dx = dy = 0.05 #Mean distribution is 40-50 micron (WR taylor et al)

#2) create a random distribution of cells and their radii
n_cells = 10
cells = zeros(n_cells, 2)
cells[:, 1] .= LinRange(0.0, 1.0, n_cells)
#cells = generate_ring_coordinates(n_cells)
#cells = rand(xmin:dx:xmax, ncells, 2)
xs = cells[:, 1]
ys = cells[:, 2]

#cells = even_map(xmin = xmin, dx = dx, xmax = xmax, ymin = ymin, dy = dy, ymax = ymax)
radii = fill(0.3, size(cells, 1)) #Switch this on to get constant radii
dist_func1(d) = ring_circle_overlap_area(d; density = 1.0, r_inner = 0.1, r_outer = 0.2, r_circle = 0.2);
cell_map = cell_map |> make_GPU

p0_dict = SAC_p0_dict()
#p0_dict["g_ACh"] = 0.215
#p0_dict["g_GABA"] = 0.0
p0 = extract_p0(p0_dict) 

u0_dict= SAC_u0_dict(mode = :PDE, ncells = n_cells)
u0 = extract_u0(u0_dict) |> CuArray{Float32}

#3) Define the problem
tspan = (0.0, 60e3)
f_PDE(du, u, p, t) = SAC_PDE(du, u, p, t, cell_map)
prob = SDEProblem(f_PDE, noise2D, u0, tspan, p0)
@time sol = solve(prob, 
     SOSRI(), #This seems to be the best solver option
     force_dtmin = true, 
     progress=true, progress_steps=1)
#save("data.jld", "initial_cond", sol[end])

#Start plotting
CUDA.allowscalar(true) #allow GPU operations to be offloaded to CPU 
Time = sol.t[1]:10:sol.t[end]
vt = hcat(map(t -> sol(t)[:,2], Time)...)|>Array
nt = hcat(map(t -> sol(t)[:,3], Time)...)|>Array
mt = hcat(map(t -> sol(t)[:,4], Time)...)|>Array
ht = hcat(map(t -> sol(t)[:,5], Time)...)|>Array
ct = hcat(map(t -> sol(t)[:,6], Time)...)|>Array
at = hcat(map(t -> sol(t)[:,7], Time)...)|>Array
bt = hcat(map(t -> sol(t)[:,8], Time)...)|>Array
et = hcat(map(t -> sol(t)[:,9], Time)...)|>Array
it = hcat(map(t -> sol(t)[:,10], Time)...)|>Array
Wt = hcat(map(t -> sol(t)[:,11], Time)...)|>Array

#Create the plot
fig1 = Figure(size = (400,800))
ax1a = Axis(fig1[1,1], title = "Voltage")
ax1b = Axis(fig1[2,1], title = "Voltage Traces")
sctV = scatter!(ax1a, cells[:,1], cells[:,2], color = sol(0.0)[:,1]|>Array, colorrange = (minimum(vt), maximum(vt)), markersize = 20.0)
for n in 1:n_cells
     lines!(ax1b, Time, vt[n, :])
end
ticker = vlines!(ax1b, [0.0])
display(fig1)

n_frames = 1000
animate_t = LinRange(0.0, sol.t[end], n_frames)
dt = animate_t[2] - animate_t[1]
fps = round(Int64, (1/dt) * 1000)

GLMakie.record(fig1, "test/SAC_model_tests/data/model_animation.mp4", animate_t, framerate = 8) do t
	println(t)
	v = sol(t)[:, 2] |> Array
	sctV.color = v
     ticker[1] = [t]
end
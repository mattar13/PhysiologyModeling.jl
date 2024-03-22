using Revise
using Pkg; Pkg.activate(".")
using ElectroPhysiology, PhysiologyModeling
Pkg.activate("test") #Activate the testing environment
using PhysiologyPlotting, GLMakie
using CUDA
CUDA.allowscalar(false)
using SparseArrays
using LinearAlgebra

#Try importing some other solver methods
import .PhysiologyModeling: CVODE_BDF, ring
import .PhysiologyModeling.DifferentialEquations: ImplicitRKMil, SKenCarp
#%%=[Run branch generation]__________________________________________________________________________________#

#1) determine the domains and spacing of cells. 
domain_x = (xmin, xmax) = (0.0, 5.0) #This is a simulation for a retina 5mm in diameter
domain_y = (ymin, ymax) = (0.0, 5.0)
dx = dy = 0.05 #Mean distribution is 40-50 micron (WR taylor et al)

#2) create a random distribution of cells and their radii
#The density of SACs in the retina is around 1200 per mm2. So if we have 5mm2 1200 * 5 = 6000
n_cells = 100 #Really pushing the model
xs = rand(xmin:dx:xmax, n_cells)
ys = rand(ymin:dy:ymax, n_cells)
connections = connect_neighbors_radius(xs, ys, 0.2) |> connection_matrix
dist_func1(d) = ring_circle_overlap_area(d; density = 0.1, r_inner = 0.1, r_outer = 0.2, r_circle = 0.2);
cell_map_CPU = CellMap(xs, ys, connections; distance_function = dist_func1);
#make sure cells are connected, if not remove unconnected cells
cell_map = cell_map_CPU |> make_GPU

# [run the model]____________________________________________________________________________#
p0_dict = SAC_p0_dict()
p0_dict["g_ACh"] = 2.0
p0_dict["g_GABA"] = 0.0
p0 = extract_p0(p0_dict) 

u0_dict= SAC_u0_dict(mode = :PDE, ncells = n_cells)
u0 = extract_u0(u0_dict) |> CuArray{Float32}

#3) Define the problem
tspan = (0.0, 60e3)
f_PDE(du, u, p, t) = SAC_PDE(du, u, p, t, cell_map)
prob = SDEProblem(f_PDE, noise2D, u0, tspan, p0)
@time sol = solve(prob, 
     #SOSRI(), #This seems to be the best solver option
     SOSRA(),
     reltol =  2e-2, abstol = 2e-2, 
     progress=true, progress_steps=1)
#save("data.jld", "initial_cond", sol[end])
println(sol)#So we can check if the solution succeeded
 
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

#%%====================================[Plot the solution]====================================#
fig1 = Figure(size = (400,800))
ax1a = Axis(fig1[1,1], title = "Calcium Imageing", xlabel = "nx", ylabel = "ny")
ax1b = Axis(fig1[2,1], title = "Calcium ROIs", xlabel = "time (ms)", ylabel = "Ct")
ax1c = Axis(fig1[3,1], title = "Voltage Traces", xlabel = "time (ms)", ylabel = "Vt (mV)")

rowsize!(fig1.layout, 1, Relative(1/2)) #Make the cell plot larger
sctV = scatter!(ax1a, cells[:,1], cells[:,2], color = sol(0.0)[:,6]|>Array, colorrange = (0.0, maximum(ct)), markersize = 20.0)
for n in 1:n_cells
     lines!(ax1b, Time, ct[n, :])
     lines!(ax1c, Time, vt[n, :])
end
tickerb = vlines!(ax1b, [0.0])
tickerc = vlines!(ax1c, [0.0])

display(fig1)

n_frames = 1000
animate_t = LinRange(0.0, sol.t[end], n_frames)
dt = animate_t[2] - animate_t[1]
fps = round(Int64, (1/dt) * 1000)

GLMakie.record(fig1, "test/SAC_model_tests/data/model_animation.mp4", animate_t, framerate = 8) do t
	println(t)
	c = sol(t)[:, 6] |> Array
	sctV.color = c
     tickerb[1] = [t]
     tickerc[1] = [t]
end
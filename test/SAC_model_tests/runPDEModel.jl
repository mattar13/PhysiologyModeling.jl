using Revise, Profile, ProfileSVG
using PhysiologyModeling
using PhysiologyPlotting
using GLMakie
using SparseArrays
import PhysiologyModeling: SRIW1, EM
using JLD

# Step 1 determine the domains and spacing of cells. 
domain_x = (xmin, xmax) = (0.0, 3.0)
domain_y = (ymin, ymax) = (0.0, 3.0)
dx = dy = 0.1 #Mean distribution is 40-50 micron (WR taylor et al)

# Step 2 create the map of cells and their radii
cells = even_map(xmin = xmin, dx = dx, xmax = xmax, ymin = ymin, dy = dy, ymax = ymax)
radii = fill(0.200, size(cells, 1)) #Switch this on to get constant radii
cell_map = CellMap(cells, radii, 
     distance_function = x -> ring_circle_overlap_area(x; density = 1.0)
);

#Step 3 create the initial conditions and parameters
#u0 = vcat(fill(vals_u0, size(cells, 1))'...) #Generate a new initial conditions
#load a previous initial conditions
u0 = load("data.jld")["initial_cond"]

#Set parameters
SAC_p0_dict["g_GABA"] = 0.0
SAC_p0_dict["I_app"] = 0.0
SAC_p0_dict["g_ACh"] = 1.0
SAC_p0_dict["g_K"] = 2.0 
p0 = extract_p0(SAC_p0_dict);

# Run the model
tspan = (0.0, 100e3)
f_PDE(du, u, p, t) = SAC_PDE(du, u, p, t, cell_map)
prob = SDEProblem(f_PDE, noise2D, u0, tspan, p0)
@time sol = solve(prob, SOSRI(), reltol = 0.1, abstol = 0.01, progress=true, progress_steps=1);

#Save the solution end. No need to warm up again
save("data.jld", "initial_cond", sol[end])

#%% Plot the figure
fDIFF = Figure(size = (1800,1000))
ax1 = Axis3(fDIFF[1,1]; aspect=(1, 1, 1))
ax2 = Axis3(fDIFF[1,2]; aspect=(1, 1, 1))
ax3 = Axis3(fDIFF[1,3]; aspect=(1, 1, 1))

surf1 = surface!(ax1, cells[:, 1], cells[:, 2], sol(0.0)[:, 1])
surf2 = surface!(ax2, cells[:, 1], cells[:, 2], sol(0.0)[:, 5])
surf3 = surface!(ax3, cells[:, 1], cells[:, 2], sol(0.0)[:, 8])

zlims!(ax1, minimum(sol[:, 1, :]), maximum(sol[:, 1, :]))
zlims!(ax2, minimum(sol[:, 5, :]), maximum(sol[:, 5, :]))
zlims!(ax3, minimum(sol[:, 8, :]), maximum(sol[:, 8, :]))

#display(fDIFF)
n_frames = 1000
animate_t = LinRange(0.0, sol.t[end], n_frames)
dt = animate_t[2] - animate_t[1]
fps = round(Int64, (1/dt) * 1000)

record(fDIFF, "test/SAC_model_tests/model_animation.mp4", animate_t, framerate = 8) do t
	println(t)
	v = sol(t)[:, 1]
     c = sol(t)[:, 5]
     e = sol(t)[:, 8]
	surf1[3] = v
     surf2[3] = c
     surf3[3] = e
     
end
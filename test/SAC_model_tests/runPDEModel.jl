using Revise
using ElectroPhysiology
using PhysiologyModeling
import PhysiologyModeling: SRIW1, EM
using Pkg; Pkg.activate("test") #Activate the testing environment
using PhysiologyPlotting, GLMakie
using CUDA
using SparseArrays
using LinearAlgebra

#%% 1) determine the domains and spacing of cells. 
domain_x = (xmin, xmax) = (0.0, 1.0)
domain_y = (ymin, ymax) = (0.0, 1.0)
dx = dy = 0.05 #Mean distribution is 40-50 micron (WR taylor et al)

#%% 2) create the map of cells and their radii
cells = even_map(xmin = xmin, dx = dx, xmax = xmax, ymin = ymin, dy = dy, ymax = ymax)
radii = fill(0.200, size(cells, 1)) #Switch this on to get constant radii
cell_map = make_GPU(CellMap(cells, radii));
u0 = vcat(fill(vals_u0, size(cells, 1))'...) |> CuArray{Float32} #Generate a new initial conditions
SAC_p0_dict["g_ACh"] = 1.0
SAC_p0_dict["g_GABA"] = 0.0
SAC_p0_dict["g_W"] = 0.075
p0 = extract_p0(SAC_p0_dict)

#3) Define the problem
tspan = (0.0, 100e3)
f_PDE(du, u, p, t) = SAC_PDE_GPU(du, u, p, t, cell_map)
prob = SDEProblem(f_PDE, noise2D, u0, tspan, p0)

#%% Pause here before running the model
@time sol = solve(prob, SOSRA(), reltol = 2e-2, abstol= 2e-2, progress=true, progress_steps=1)
sol.t

#Save the solution end. No need to warm up again
#save("data.jld", "initial_cond", sol[end])
#%% Find the zlims
zlims1 =  (minimum(sol[:, 1, :]), maximum(sol[:, 1, :]))
zlims2 =  (minimum(sol[:, 8, :]), maximum(sol[:, 8, :]))
zlims3 =  (minimum(sol[:, 9, :]), maximum(sol[:, 9, :]))
zlims4 =  (minimum(sol[:, 5, :]), maximum(sol[:, 5, :]))
zlims5 =  (minimum(sol[:, 6, :]), maximum(sol[:, 6, :]))
zlims6 =  (minimum(sol[:, 7, :]), maximum(sol[:, 7, :]))
#%% Plot the figure
fDIFF = Figure(size = (1800,1000))
ax1 = Axis3(fDIFF[1,1]; aspect=(1, 1, 1), title = "Voltage")
ax2 = Axis3(fDIFF[1,2]; aspect=(1, 1, 1), title = "Acetylcholine")
ax3 = Axis3(fDIFF[1,3]; aspect=(1, 1, 1), title = "GABA")

ax4 = Axis3(fDIFF[2,1]; aspect=(1, 1, 1), title = "Calcium")
ax5 = Axis3(fDIFF[2,2]; aspect=(1, 1, 1), title = "cAMP")
ax6 = Axis3(fDIFF[2,3]; aspect=(1, 1, 1), title = "TREK1")

surf1 = surface!(ax1, cells[:, 1], cells[:, 2], sol(0.0)[:, 1])
surf2 = surface!(ax2, cells[:, 1], cells[:, 2], sol(0.0)[:, 8])
surf3 = surface!(ax3, cells[:, 1], cells[:, 2], sol(0.0)[:, 9])

surf4 = surface!(ax4, cells[:, 1], cells[:, 2], sol(0.0)[:, 5])
surf5 = surface!(ax5, cells[:, 1], cells[:, 2], sol(0.0)[:, 6])
surf6 = surface!(ax6, cells[:, 1], cells[:, 2], sol(0.0)[:, 7])

zlims!(ax1, zlims1)
zlims!(ax2, zlims2)
zlims!(ax3, zlims3)
zlims!(ax4, zlims4)
zlims!(ax5, zlims5)
zlims!(ax6, zlims6)

#display(fDIFF)
n_frames = 1000
animate_t = LinRange(0.0, sol.t[end], n_frames)
dt = animate_t[2] - animate_t[1]
fps = round(Int64, (1/dt) * 1000)

GLMakie.record(fDIFF, "test/SAC_model_tests/model_animation.mp4", animate_t, framerate = 8) do t
	println(t)
	v = sol(t)[:, 1]
     e = sol(t)[:, 8]
     i = sol(t)[:, 9]
     c = sol(t)[:, 5]
     a = sol(t)[:, 6]
     b = sol(t)[:, 7]
	surf1[3] = v
     surf2[3] = e
     surf3[3] = i  
     surf4[3] = c
     surf5[3] = a
     surf6[3] = b
end
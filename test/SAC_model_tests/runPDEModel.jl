using Revise
using ElectroPhysiology
using PhysiologyModeling
import PhysiologyModeling: SRIW1, EM, ring
using Pkg; Pkg.activate("test") #Activate the testing environment
using PhysiologyPlotting, GLMakie
using CUDA
CUDA.allowscalar(false)
using SparseArrays
using LinearAlgebra

#%% 1) determine the domains and spacing of cells. 
domain_x = (xmin, xmax) = (0.0, 1.0)
domain_y = (ymin, ymax) = (0.0, 1.0)
dx = dy = 0.05 #Mean distribution is 40-50 micron (WR taylor et al)

#2) create the map of cells and their radii
cells = even_map(xmin = xmin, dx = dx, xmax = xmax, ymin = ymin, dy = dy, ymax = ymax)
radii = fill(0.2, size(cells, 1)) #Switch this on to get constant radii
dist_func(d) = ring(d; max_strength = 0.005, max_dist = 0.15)
cell_map = CellMap(cells, radii; distance_function = dist_func) |> make_GPU;

ics = extract_u0(SAC_u0_dict)
u0_CPU = vcat(fill(ics, size(cells, 1))'...) 
u0_CPU[1, 1] = 0.0
u0 = u0_CPU |> CuArray{Float32} #Generate a new initial conditions

SAC_p0_dict["g_ACh"] = 1.0
SAC_p0_dict["g_GABA"] = 0.0
SAC_p0_dict["g_W"] = 0.075
p0 = extract_p0(SAC_p0_dict) 

#3) Define the problem
tspan = (0.0, 100e3)
f_PDE(du, u, p, t) = SAC_PDE_GPU(du, u, p, t, cell_map)
prob = SDEProblem(f_PDE, noise2D, u0, tspan, p0)

# Pause here before running the model
@time sol = solve(prob, SOSRI(), reltol = 2e-1, abstol= 2e-1, progress=true, progress_steps=1)
#save("data.jld", "initial_cond", sol[end])
sol.t

#%% Find the zlims
CUDA.allowscalar(true) #
zlims1 =  (minimum(sol[:, 1, :]), maximum(sol[:, 1, :]))
zlims2 =  (minimum(sol[:, 8, :]), maximum(sol[:, 8, :]))
zlims3 =  (minimum(sol[:, 9, :]), maximum(sol[:, 9, :]))
zlims4 =  (minimum(sol[:, 5, :]), maximum(sol[:, 5, :]))
zlims5 =  (minimum(sol[:, 6, :]), maximum(sol[:, 6, :]))
zlims6 =  (minimum(sol[:, 7, :]), maximum(sol[:, 7, :]))

# Plot the figure
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

#%%
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
for i in rand(1:size(sol, 1), 10)
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
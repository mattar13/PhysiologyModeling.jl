using Revise, Profile, ProfileSVG
using PhysiologyModeling
using PhysiologyPlotting
using GLMakie
using SparseArrays
import PhysiologyModeling: SRIW1, EM

# Step 1 determine the domains and spacing of cells. 
domain_x = (xmin, xmax) = (0.0, 1.0)
domain_y = (ymin, ymax) = (0.0, 1.0)
dx = dy = 0.057 #Mean distribution is 40-50 micron (WR taylor et al)

# Step 2 create the map of cells and their radii
cells = even_map(xmin = xmin, dx = dx, xmax = xmax, ymin = ymin, dy = dy, ymax = ymax)
radii = fill(0.200, size(cells, 1)) #Switch this on to get constant radii
cell_map = CellMap(cells, radii);

#Step 3 create the initial conditions and parameters
u0 = vcat(fill(vals_u0, size(cells, 1))'...)

#Set parameters
SAC_p0_dict["g_GABA"] = 0.0
SAC_p0_dict["I_app"] = 0.0
SAC_p0_dict["g_ACh"] = 1.0
p0 = extract_p0(SAC_p0_dict);

# Run the model
tspan = (0.0, 60.0)
f_PDE(du, u, p, t) = SAC_PDE(du, u, p, t, cell_map)
prob = SDEProblem(f_PDE, noise2D, u0, tspan, p0)
@profile sol = solve(prob, SOSRI(), reltol = 0.01, abstol = 0.01, progress=true, progress_steps=1);

ProfileSVG.save("profile.svg")

#%% Plot the figure
fDIFF = Figure(resolution = (1800,1000))
ax1 = Axis3(fDIFF[1,1]; aspect=(1, 1, 1))
ax2 = Axis3(fDIFF[1,2]; aspect=(1, 1, 1))
ax3 = Axis3(fDIFF[1,3]; aspect=(1, 1, 1))
ax4 = Axis3(fDIFF[2,1]; aspect=(1, 1, 1))
ax5 = Axis3(fDIFF[2,2]; aspect=(1, 1, 1))
ax6 = Axis3(fDIFF[2,3]; aspect=(1, 1, 1))

surf1 = surface!(ax1, cells[:, 1], cells[:, 2], sol(100.0)[:, 1])
surf2 = surface!(ax2, cells[:, 1], cells[:, 2], sol(100.0)[:, 10])
surf3 = surface!(ax3, cells[:, 1], cells[:, 2], sol(100.0)[:, 8])
surf4 = surface!(ax4, cells[:, 1], cells[:, 2], sol(100.0)[:, 5])
surf5 = surface!(ax5, cells[:, 1], cells[:, 2], sol(100.0)[:, 6])
surf6 = surface!(ax6, cells[:, 1], cells[:, 2], sol(100.0)[:, 7])

zlims!(ax1, minimum(sol[:, 1, :]), maximum(sol[:, 1, :]))
zlims!(ax2, minimum(sol[:, 10, :]), maximum(sol[:, 10, :]))
zlims!(ax3, minimum(sol[:, 8, :]), maximum(sol[:, 8, :]))
zlims!(ax4, minimum(sol[:, 5, :]), maximum(sol[:, 5, :]))
zlims!(ax5, minimum(sol[:, 6, :]), maximum(sol[:, 6, :]))
zlims!(ax6, minimum(sol[:, 7, :]), maximum(sol[:, 7, :]))
#zlims!(ax7, minimum(sol[:, 7, :]), maximum(sol[:, 7, :]))
#zlims!(ax8, minimum(sol[:, 8, :]), maximum(sol[:, 8, :]))
#zlims!(ax9, minimum(sol[:, 9, :]), maximum(sol[:, 9, :]))
#zlims!(ax10, minimum(sol[:, 10, :]), maximum(sol[:, 10, :]))

#display(fDIFF)
#
n_frames = 1000
animate_t = LinRange(0.0, sol.t[end], n_frames)
dt = animate_t[2] - animate_t[1]
fps = round(Int64, (1/dt) * 1000)

record(fDIFF, "test/SAC_model_tests/model_animation.mp4", animate_t, framerate = 8) do t
	println(t)
	v = sol(t)[:, 1]
     n = sol(t)[:, 2]
     m = sol(t)[:, 3]
     h = sol(t)[:, 4]
     c = sol(t)[:, 5]
     a = sol(t)[:, 6]
     b = sol(t)[:, 7]
     e = sol(t)[:, 8]
     i = sol(t)[:, 9]
     W = sol(t)[:, 10]
	surf1[3] = v
     surf2[3] = W
     surf3[3] = e
     surf4[3] = c
     surf5[3] = a
     surf6[3] = b
     #surf7[3] = b
     #surf8[3] = e
     #surf9[3] = i
     #surf10[3] = W
end
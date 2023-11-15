using Revise
using PhysiologyModeling
using PhysiologyPlotting
using GLMakie
using SparseArrays

# Step 1 determine the domains and spacing of cells. 
domain_x = (xmin, xmax) = (0.0, 2.0)
domain_y = (ymin, ymax) = (0.0, 2.0)
dx = dy = 0.05 #Mean distribution is 40-50 micron (WR taylor et al)

# Step 2 create the map of cells and their radii
#cells = rand(50, 2) .|> abs
#radii = rand(0.05:0.01:0.180, size(cells, 1)) #Switch this on to get random radii
cells = even_map(xmin = xmin, dx = dx, xmax = xmax, ymin = ymin, dy = dy, ymax = ymax)
radii = fill(0.18, size(cells, 1)) #Switch this on to get constant radii
cell_map = CellMap(cells, radii, max_strength = 0.005);

#Step 3 create the initial conditions and parameters
u0 = vcat(fill(vals_u0, size(cells, 1))'...)
mid = round(Int64, size(cell_map.xs, 1)/2)+1
#u0[mid, 8] = 5.0

#Set parameters
SAC_p0_dict["g_GABA"] = 0.0
SAC_p0_dict["I_app"] = 0.0
SAC_p0_dict["g_ACh"] = 1.0
MAP_p = (cell_map, extract_p0(SAC_p0_dict));

#%% Warmup the problem for 50s
probWARM = SDEProblem(SAC_PDE, noise2D, u0, (0.0, 60e3), MAP_p)
@time solWARM = solve(probWARM, SOSRI(), save_everystep = false, save_start = false, reltol=0.1, abstol=0.1, progress=true, progress_steps=1)
solWARM[end]

#%% Run the model
tspan = (0.0, 60e3)
#prob = ODEProblem(SAC_PDE, u0, tspan, MAP_p)
prob = SDEProblem(SAC_PDE, noise2D, solWARM[end], tspan, MAP_p)
#sol = solve(prob, TRBDF2(), reltol=0.01, abstol=0.01, progress=true, progress_steps=1)
@time sol = solve(prob, SOSRI(), saveat = 1.0, reltol=0.1, abstol=0.1, progress=true, progress_steps=1)

#%% Plot the figure
fDIFF = Figure(resolution = (1800,1000))
ax1 = Axis3(fDIFF[1,1]; aspect=(1, 1, 1))
ax2 = Axis3(fDIFF[1,2]; aspect=(1, 1, 1))
ax3 = Axis3(fDIFF[1,3]; aspect=(1, 1, 1))
ax4 = Axis3(fDIFF[2,1]; aspect=(1, 1, 1))
ax5 = Axis3(fDIFF[2,2]; aspect=(1, 1, 1))
ax6 = Axis3(fDIFF[2,3]; aspect=(1, 1, 1))
#ax7 = Axis3(fDIFF[2,2]; aspect=(1, 1, 1))
#ax8 = Axis3(fDIFF[2,3]; aspect=(1, 1, 1))
#ax9 = Axis3(fDIFF[2,4]; aspect=(1, 1, 1))
#ax10 = Axis3(fDIFF[2,5]; aspect=(1, 1, 1))

surf1 = surface!(ax1, cells[:, 1], cells[:, 2], sol(100.0)[:, 1])#sol(0.0)[:, 1])
surf2 = surface!(ax2, cells[:, 1], cells[:, 2], sol(100.0)[:, 10])
surf3 = surface!(ax3, cells[:, 1], cells[:, 2], sol(100.0)[:, 8])
surf4 = surface!(ax4, cells[:, 1], cells[:, 2], sol(100.0)[:, 5])
surf5 = surface!(ax5, cells[:, 1], cells[:, 2], sol(100.0)[:, 6])
surf6 = surface!(ax6, cells[:, 1], cells[:, 2], sol(100.0)[:, 7])
#surf7 = surface!(ax7, cells[:, 1], cells[:, 2], sol(100.0)[:, 7])
#surf8 = surface!(ax8, cells[:, 1], cells[:, 2], sol(100.0)[:, 8])
#surf9 = surface!(ax9, cells[:, 1], cells[:, 2], sol(100.0)[:, 9])
#surf10 = surface!(ax10, cells[:, 1], cells[:, 2], sol(100.0)[:, 10])

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


sol_arr[:, 8, :] |> maximum

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
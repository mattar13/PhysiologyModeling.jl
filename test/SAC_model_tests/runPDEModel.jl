using Revise
using PhysiologyModeling
using PhysiologyPlotting
using GLMakie
using SparseArrays

#%% Step 1 determine the domains and spacing of cells. 
domain_x = (xmin, xmax) = (0.0, 1.0)
domain_y = (ymin, ymax) = (0.0, 1.0)
dx = dy = 0.05 #Mean distribution is 40-50 micron (WR taylor et al)

# Step 2 create the map of cells and their radii
cells = even_map(xmin = xmin, dx = dx, xmax = xmax, ymin = ymin, dy = dy, ymax = ymax)
#radii = rand(0.05:0.01:0.180, size(cells, 1)) #Switch this on to get random radii
radii = fill(0.180, size(cells, 1)) #Switch this on to get constant radii
cell_map = CellMap(cells, radii, strength = 0.001);

#Step 3 create the initial conditions and parameters
u0 = vcat(fill(vals_u0, size(cells, 1))'...)
mid = round(Int64, size(cell_map.xs, 1)/2)+1
u0[mid, 8] = 5.0

#Set parameters
SAC_p0_dict["g_GABA"] = 0.0
SAC_p0_dict["I_app"] = 0.0
MAP_p = (cell_map, extract_p0(SAC_p0_dict));

tspan = (0.0, 1000.0)
prob = ODEProblem(SAC_PDE, u0, tspan, MAP_p)
sol = solve(prob, progress=true, progress_steps=1)

#%% Plot the voltage, calcium and ACh
time = sol.t
fODE = Figure(resolution = (800, 800))
ax1 = Axis(fODE[1,1])
ax2 = Axis(fODE[2,1])
ax3 = Axis(fODE[3,1])
ax4 = Axis(fODE[4,1])
for i in axes(u0, 1)
     println(i)
     lines!(ax1, time, map(t -> sol(t)[i, 1], time))
     lines!(ax2, time, map(t -> sol(t)[i, 5], time))
     lines!(ax3, time, map(t -> sol(t)[i, 8], time))
     lines!(ax4, time, map(t -> sol(t)[i, 9], time))
end
display(fODE)

#%%
fDIFF = Figure(resolution = (1000,1000))



#%% Plot the figure
ax11 = Axis3(fDIFF[1,1]; aspect=(1, 1, 1))
ax12 = Axis3(fDIFF[1,2]; aspect=(1, 1, 1))

ax21 = Axis(fDIFF[2,1])
ax22 = Axis(fDIFF[2,2])

surface!(ax11, cell_map.xs, cell_map.ys, sol(0.0))
scatter!(ax21, cell_map.xs, cell_map.ys, color = sol(0.0), markersize = 15.0)
scatter!(ax21, cell_map.xs[mid], cell_map.ys[mid], strokewidth = 1, color = :transparent, markersize = cell_map.radius[mid]*2, markerspace = :data)
surf2 = surface!(ax12, cells[:, 1], cells[:, 2], sol(0.0), colorrange = (0.0, 1.0), markersize = 35.0)
scat2 = scatter!(ax22, cells[:, 1], cells[:, 2], color = sol(0.0), colorrange = (0.0, 1.0), markersize = 15.0)
scatter!(ax22, cell_map.xs[mid], cell_map.ys[mid], strokewidth = 1, color = :transparent, markersize = cell_map.radius[mid]*2, markerspace = :data)
xlims!(ax11, (xmin, xmax))
ylims!(ax11, (ymin, ymax))
zlims!(ax11, (0.0, maximum(u0)))

xlims!(ax12, (xmin, xmax))
ylims!(ax12, (ymin, ymax))
zlims!(ax12, (0.0, maximum(u0)))
record(fDIFF, "time_animation.mp4", 0.0:10.0:sol.t[end]) do t
	println(t)
	u = sol(t)
	surf2[3] = u
	scat2.color = u
end
using Revise
using PhysiologyModeling
using PhysiologyPlotting
using GLMakie

#%% Set up plotting

# Lets try some advanced mapping
n_cells = 100
n_epoch = 1000
k = 10 #Each cell will have 5 nearest neighbors
min_vs = zeros(n_epoch+1)
max_vs = zeros(n_epoch+1)
sum_vs = zeros(n_epoch+1)
avg_vs = zeros(n_epoch+1)
#raster = 

cell_map = CellMap(n_cells, k; strength = 0.005, map_spacing = :even, dx = 0.01, dy = 0.01);
u = zeros(n_cells)
val_record = zeros(n_cells, n_epoch+1)

u[5] = 100.0
u[25] = 100.0
u[95] = 100.0


#Change plotting to Makie
#xmap, ymap, grid_pre = rasterize(cell_map, dx = 0.02, dy = 0.02)

min_vs[1] = minimum(u)
max_vs[1] = maximum(u)
sum_vs[1] = sum(u)
avg_vs[1] = sum(u)/n_cells
val_record[:, 1] .= u

for i in 1:n_epoch
	println(i)
	du = zeros(size(u))
	#try adding some iteration in
	∇α(du, u, cell_map)
	u .+= du

	min_vs[i+1] = minimum(u)
	max_vs[i+1] = maximum(u)
	sum_vs[i+1] = sum(u)
	avg_vs[i+1] = sum(u)/n_cells
	val_record[:, i+1] .= u 
end
#xmap, ymap, grid_post = rasterize(cell_map)

#%%
f = Figure(resolution =(1600, 1600))
ax11 = Axis(f[1,1])
ax12 = Axis(f[1,2])
ax21 = Axis(f[2,1])
ax22 = Axis(f[2,2])
ax31 = Axis(f[3,1])
ax32 = Axis(f[3,2])

lines!(ax21, 1:n_epoch+1, sum_vs)
lines!(ax22, 1:n_epoch+1, avg_vs)
lines!(ax22, 1:n_epoch+1, min_vs)
lines!(ax22, 1:n_epoch+1, max_vs)
scatter!(ax11, cell_map.xs, cell_map.ys, color = val_record[:, 1], markersize = 35.0)
scat_plot = scatter!(ax12, cell_map.xs, cell_map.ys, color = val_record[:, 1], markersize = 35.0)
#heatmap!(ax31, xmap, ymap, grid_pre)
#hmp = heatmap!(ax32, xmap, ymap, grid_post)

record(f, "animation.mp4", 1:n_epoch+1, framerate = 5) do t
	#println(t)
	scat_plot.color = val_record[:, t]
	#hmp.values = 
end
using Revise
using Pkg; Pkg.activate(".")
using PhysiologyModeling
import PhysiologyModeling: ring, Φe, IACh, IGABA, ħe, ħi, ring_circle_overlap_area
Pkg.activate("test")
using PhysiologyPlotting
using GLMakie
using SparseArrays
import .PhysiologyModeling.euclidean_distance

#%% [Draw the maps for distance heatmaps] ____________________________________________#
origin = (0.0, 0.0)
xrng = LinRange(-0.750, 0.750, 100)
yrng = LinRange(-0.750, 0.750, 100)
coords = Iterators.product(xrng, yrng) |> collect
distances = zeros(size(coords)...)
for ix in axes(coords, 1), iy in axes(coords,2)
     distances[ix, iy] = euclidean_distance(origin, coords[ix, iy])
end
dist_func1(d) = ring_circle_overlap_area(d; density = 0.01, r_inner = 0.1, r_outer = 0.18, r_circle = 0.18)
dist_func2(d) = ring(d; density = 0.01, max_dist = 0.18, slope = 0.025)
dist_func3(d) = d
strengths1 = dist_func1.(distances)
strengths2 = dist_func2.(distances)
strengths3 = dist_func3.(distances)
distances

# [Plot the function] _______________________________________________________________#
fig1 = Figure(size = (1200, 300))
ax1a = Axis(fig1[1,1])
ax1b = Axis(fig1[1,2])
ax1c = Axis(fig1[1,3])
ax1d = Axis(fig1[1,4])
heatmap!(ax1a, xrng, yrng, distances)
heatmap!(ax1b, xrng, yrng, strengths1)
heatmap!(ax1c, xrng, yrng, strengths2)
heatmap!(ax1d, xrng, yrng, strengths3)
scatter!(ax1a, origin); scatter!(ax1b, origin); scatter!(ax1c, origin)
save("test/SAC_model_tests/data/MappingFunc.png", fig1)

#%% [Drawing random cell maps] ________________________________________________________#
#1) determine the domains and spacing of cells. 
domain_x = (xmin, xmax) = (0.0, 5.0) #This is a simulation for a retina 5mm in diameter
domain_y = (ymin, ymax) = (0.0, 5.0)
dx = dy = 0.05 #Mean distribution is 40-50 micron (WR taylor et al)

#2) create a random distribution of cells and their radii
#The density of SACs in the retina is around 1200 per mm2. So if we have 5mm2 1200 * 5 = 6000
n_cells = 600 #Really pushing the model
xs, ys = create_random_map(n_cells, dx = dy = 0.01, xmax = 2.0, ymax = 2.0)

connection_list = connect_neighbors_radius(xs, ys, 0.2)
connection_list = connection_list[findall(c -> c[3] > 0.0, connection_list)]
connections = connection_matrix(connection_list, m = length(xs), n = length(ys))
dist_func1(d) = ring(d; density = 0.01, max_dist = 0.18, slope = 0.025);
cell_map_CPU = CellMap(xs, ys, connections; distance_function = dist_func1);
cell_map_CPU.strength[2,:]

fig2 = Figure(size = (400,400))
ax2a = Axis(fig2[1,1], title = "Coordinates", xlabel = "nx", ylabel = "ny")
#ax2b = Axis(fig2[2,1], title = "Calcium ROIs", xlabel = "time (ms)", ylabel = "Ct")
#ax2c = Axis(fig2[3,1], title = "Voltage Traces", xlabel = "time (ms)", ylabel = "Vt (mV)")
rows, cols, data = findnz(cell_map_CPU.strength)
for i in eachindex(rows)
    r = rows[i]
    c = cols[i]
    d = data[i]
    lines!(ax2a, [xs[r], xs[c]], [ys[r], ys[c]], colormap = :viridis, color = [d], colorrange = (minimum(data), maximum(data)))
end
sctV = scatter!(ax2a, xs, ys, markersize = 15.0)
display(fig2)
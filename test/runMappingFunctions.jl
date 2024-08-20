using ElectroPhysiology, PhysiologyModeling
using Pkg; Pkg.activate("test")
using PhysiologyPlotting, GLMakie
import PhysiologyModeling: euclidean_distance, calculate_linear_bias, find_angle

#%% [Eventually we want to generalize the mapping function] ____________________________________________#
origin = (0.0, 0.0)
xs = LinRange(-100.0, 100.0, 1000)
ys = LinRange(-100.0, 100.0, 1000)
coords = Iterators.product(xs, ys) |> collect

rinc_str = zeros(size(coords)...)
ring_str = zeros(size(coords)...)
cons_str = zeros(size(coords)...)
angles = zeros(size(coords)...)
dist_func1(p1, p2) = RING_CIRC(p1, p2; density = 80e-6, r_inner = 40.0, r_outer = 80.0, r_circle = 80.0) 
dist_func2(p1, p2) = RING_CIRC(p1, p2; density = 80e-6, r_inner = 40.0, r_outer = 80.0, r_circle = 80.0)

angle_bias = 180
for ix in axes(coords, 1), iy in axes(coords,2)
    bias = calculate_exponential_bias(find_angle(origin, coords[ix, iy]), 90.0, decay_rate = 0.03)
    rinc_str[ix, iy] = dist_func1(origin, coords[ix, iy])
    ring_str[ix, iy] = dist_func2(origin, coords[ix, iy]) * bias
end

# [Plot the function] _______________________________________________________________#
fig1 = Figure(size = (600, 500))
ax1a = Axis(fig1[1,1], title = "RingCirc")
ax1b = Axis(fig1[1,2], title = "RingCircBias")

ax2a = Axis(fig1[3,1], xlabel = "Vert. Distance (mm)", ylabel = "Strength")
ax2b = Axis(fig1[3,2], xlabel = "Vert. Distance (mm)", ylabel = "Strength")

hma = heatmap!(ax1a, xs, ys, rinc_str)
hmb = heatmap!(ax1b, xs, ys, ring_str)
Colorbar(fig1[2, 1], hma, vertical = false, label = "Fluorescence (px)")
Colorbar(fig1[2, 2], hmb, vertical = false, label = "Fluorescence (px)")

lines!(ax2a, xs, sum(rinc_str, dims = 2)[:,1], color = :black)
lines!(ax2b, xs, sum(ring_str, dims = 2)[:,1], color = :black)

hma = heatmap!(ax1a, xs, ys, rinc_str)
hmb = heatmap!(ax1b, xs, ys, ring_str)
Colorbar(fig1[2, 1], hma, vertical = false, label = "Fluorescence (px)")
Colorbar(fig1[2, 2], hmb, vertical = false, label = "Fluorescence (px)")

scatter!(ax1a, origin); scatter!(ax1b, origin)
rowsize!(fig1.layout, 2, 0.12)
rowsize!(fig1.layout, 3, 100.00)
display(fig1)

#save("Modeling/Results/MappingFunc.png", fig1)

#%% [Changing how coordinates are determined]
#1) determine the domains and spacing of cells. 
domain_x = (xmin, xmax) = (-5.0, 5.0) #This is a simulation for a retina 5mm in diameter
domain_y = (ymin, ymax) = (-5.0, 5.0)
dx = dy = 0.05 #Mean distribution is 40-50 micron (WR taylor et al)

#2) create a random distribution of cells and their radii
#The density of SACs in the retina is around 1200 per mm2. So if we have 5mm2 1200 * 5 = 6000
n_cells = 600 #Really pushing the model
xs, ys = create_random_map(n_cells, 
     xmin = xmin, dx = dx, xmax = xmax, 
     ymin = ymin, dy = dy, ymax = ymax
)
connection_list = connect_neighbors_radius(xs, ys, 0.18, self_connecting = false)
connections = connection_matrix(connection_list, m = length(xs), n = length(ys))

#Plot your bias functions
bias_func(p1, p2) = calculate_exponential_bias(find_angle(p1, p2), 90.0)
dist_func(p1, p2) = RING_CIRC(p1, p2) * bias_func(p1, p2)
cells = CellMap(xs, ys, connections, distance_function = dist_func);

#Reset the first point at each spot 
fig2 = Figure(size = (800, 800))
ax1 = Axis(fig2[1,1], 
    xlabel = "X distance", ylabel = "Y distance"
)
scatter!(ax1, (0.0, 0.0), color = :blue, markersize = 20.0)

rows, cols, data = findnz(cells.strength)
for i in eachindex(rows)
    r = rows[i]
    c = cols[i]
    d = data[i]
    p1 = (xs[r], ys[r])# .- p2
    p2 = (xs[c], ys[c]) .- p1
    scatter!(ax1, p2, color = [d], colorrange = (0.0, maximum(values)), colormap = :viridis, markersize = 50.0)
end

display(fig2)
save("test/Modeling/Results/ScatterMappingFunc.png", fig2)
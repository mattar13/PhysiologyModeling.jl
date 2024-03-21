using Revise
using ElectroPhysiology, PhysiologyModeling
using Pkg; Pkg.activate("test") #Activate the testing environment
using PhysiologyPlotting, GLMakie
using SparseArrays

#%%=[Run branch generation]__________________________________________________________________________________#

xs, ys, connection_list = create_dendrogram_map(4, 2, 5, radius = 0.01, branch_distance = 0.65)
connection_list
xs .+= rand(length(xs))/1000
xs .+= rand(length(ys))/100
connections = connection_matrix(connection_list)

dist_func1(d) = ring_circle_overlap_area(d; density = 0.1, r_inner = 0.1, r_outer = 0.2, r_circle = 0.2);
cell_map_CPU = CellMap(xs, ys, connections; distance_function = dist_func1);
#make sure cells are connected, if not remove unconnected cells
cell_map = cell_map_CPU |> make_GPU
cell_map_CPU.strength[:, 2]

for connection in connection_list
    println(connection)
end

#%% [Plot the solution]____________________________________________________________________________________________________________#
rows, cols = findnz(connections)

fig = Figure()
ax1 = Axis(fig[1,1])
scatter!(ax1, xs, ys, colormap = :viridis, color = 1:length(xs))
for (r, c) in zip(rows, cols)
    #println(r)
    #println(c)
    lines!(ax1, [xs[r], xs[c]], [ys[r], ys[c]], color = :black)
end
display(fig)
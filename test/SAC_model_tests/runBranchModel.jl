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
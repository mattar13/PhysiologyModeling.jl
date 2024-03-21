using Revise
using ElectroPhysiology, PhysiologyModeling
using Pkg; Pkg.activate("test") #Activate the testing environment
using PhysiologyPlotting, GLMakie
using SparseArrays

#%%=[Run branch generation]__________________________________________________________________________________#

xs, ys, conns = create_dendrogram_map(4, 2, 5, radius = 0.01, branch_distance = 0.65)
xs .+= rand(length(xs))/1000
xs .+= rand(length(ys))/1000
connection_matrix = create_connection_matrix(conns)

#%% [Plot the solution]____________________________________________________________________________________________________________#
rows, cols = findnz(connection_matrix)

fig = Figure()
ax1 = Axis(fig[1,1])
scatter!(ax1, xs, ys, colormap = :viridis, color = 1:length(xs))
for (r, c) in zip(rows, cols)
    #println(r)
    #println(c)
    lines!(ax1, [xs[r], xs[c]], [ys[r], ys[c]], color = :black)
end
display(fig)

a = [1,2,3,4]
repeat(a, inner = 2)
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

strengths1 = dist_func1.(distances)
strengths2 = dist_func2.(distances)
distances

#%% [Drawing random cell maps] ________________________________________________________#
xs, ys = create_random_map(300)


#%% [Plot the function] _______________________________________________________________#
fmap = Figure(size = (900, 600))
axA1 = Axis(fmap[1,1])
axA2 = Axis(fmap[1,2])
axA3 = Axis(fmap[1,3])
heatmap!(axA1, xrng, yrng, distances)
heatmap!(axA2, xrng, yrng, strengths1)
heatmap!(axA3, xrng, yrng, strengths2)
scatter!(axA1, origin); scatter!(axA2, origin); scatter!(axA3, origin)

axB1 = Axis(fmap[2,1])
axB2 = Axis(fmap[2,2])
axB3 = Axis(fmap[2,3])

display(fmap)
#%%
save("test/SAC_model_tests/data/DistanceFunc.png", f_DIST)
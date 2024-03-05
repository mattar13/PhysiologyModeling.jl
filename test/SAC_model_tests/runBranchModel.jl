using Revise
using ElectroPhysiology, PhysiologyModeling
using Pkg; Pkg.activate("test") #Activate the testing environment
using PhysiologyPlotting, GLMakie
#This uses all of the same auxillary equations

import .PhysiologyModeling.generate_circles

n = 4
z = 2
l = 2
radius_increment = 0.005
circ = generate_circles(n, z, l)
circ = hcat(circ...)'
#Try this for radial lines


fig = Figure(size = (800,800))
ax1 = Axis(fig[1,1])
scatter!(ax1, circ)
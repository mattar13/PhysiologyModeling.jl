using Revise, Profile, ProfileSVG
using PhysiologyModeling
using PhysiologyPlotting
using GLMakie
using SparseArrays, LinearAlgebra
import PhysiologyModeling: Φe

#1) determine the domains and spacing of cells. 
domain_x = (xmin, xmax) = (0.0, 0.1)
domain_y = (ymin, ymax) = (0.0, 0.1)
dx = dy = 0.05 #Mean distribution is 40-50 micron (WR taylor et al)

#2) create the map of cells and their radii
cells = even_map(xmin = xmin, dx = dx, xmax = xmax, ymin = ymin, dy = dy, ymax = ymax)
radii = fill(0.200, size(cells, 1)) #Switch this on to get constant radii
cell_map = CellMap(cells, radii);

#3) Define the initial state and timespan
u0 = zeros(size(cell_map.connections, 1))
mid = round(Int64, size(cell_map.connections, 1)/2)+1
u0[mid] = 1.0

#4) Run model
tspan = (0.0, 5000)
f_diffuse(du, u, p, t) = DIFFUSION_MODEL(du, u, p, t; active_cell = mid)
probSDE = SDEProblem(f_diffuse, DIFFUSION_NOISE, u0, tspan, cell_map)
@time begin
     @profile sol = solve(probSDE, SOSRI(), reltol=0.01, abstol=0.01, progress=true, progress_steps=1)
end
ProfileSVG.save("diffusion_profile.svg")

#%% Better calculations
u = u0
du = similar(u)
flow_in = similar(du)
flow_out = similar(du)
mul!(flow_in, cell_map.strength, u)
flow_out = cell_map.strength_out .* u 
du = flow_out + flow_in
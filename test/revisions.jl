using Revise
using ElectroPhysiology, PhysiologyModeling
using Pkg;Pkg.activate("test")
using DifferentialEquations
using GLMakie

#%% 1) determine the domains and spacing of cells. 
domain_x = (xmin, xmax) = (0.0, 0.190) #This is a simulation for a retina 5mm in diameter
domain_y = (ymin, ymax) = (0.0, 0.190)
dx = dy = 0.005 #Mean distribution is 40-50 micron (WR taylor et al)

#2) create a random distribution of cells and their radii
#The density of SACs in the retina is around 1200 per mm2. So if we have 5mm2 1200 * 5 = 6000
n_cells = 200 #Really pushing the model
xs, ys = create_random_map(n_cells, 
     xmin = xmin, dx = dx, xmax = xmax, 
     ymin = ymin, dy = dy, ymax = ymax
)
connection_list = connect_neighbors_radius(xs, ys, 0.18, self_connecting = false)
connections = connection_matrix(connection_list, m = length(xs), n = length(ys))

ACH_dist_func(p1, p2) = RING_CIRC(p1, p2) 
cell_map_ACH = CellMap(xs, ys, connections; distance_function = ACH_dist_func) |> make_GPU;
#Make the map for GABA
bias_func(p1, p2) = calculate_exponential_bias(find_angle(p1, p2), 90.0)
GABA_dist_func(p1, p2) = RING_CIRC(p1, p2; density = 0.05) * bias_func(p1, p2)
cell_map_GABA = CellMap(xs, ys, connections; distance_function = GABA_dist_func) |> make_GPU;

#[run the model]____________________________________________________________________________#
#Load parameters
p0_dict = SAC_p0_dict()
p0 = extract_p0(p0_dict) 

#Load initial conditions
u0_dict= SAC_u0_dict(mode = :PDE, n_cells = n_cells)
u0 = extract_u0(u0_dict) |> CuArray{Float32}

#3) Define the problem
tspan = (0.0, 300e3)
dt = 1.0
tseries = tspan[1]:dt:tspan[2]

f_PDE(du, u, p, t) = SAC_PDE(du, u, p, t, cell_map_ACH, cell_map_GABA) #for now diffusion is the same in both directions
prob = SDEProblem(f_PDE, noise2D, u0, tspan, p0)
@time sol = solve(prob, 
     #SOSRI(), #This seems to be the best solver option
     SOSRA(),
     #tstops = tseries, callback = cb,
     reltol =  0.1, abstol = 0.1, 
     progress=true, progress_steps=1
)
 

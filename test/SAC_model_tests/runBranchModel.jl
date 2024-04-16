using Revise
using Pkg; Pkg.activate(".")
using ElectroPhysiology, PhysiologyModeling
Pkg.activate("test") #Activate the testing environment
using PhysiologyPlotting, GLMakie
using CUDA
CUDA.allowscalar(false)
using SparseArrays
using LinearAlgebra

save_fn = ""

#%%=[Run branch generation]__________________________________________________________________________________#
radial = 4
branches = 2
layers = 3

xs, ys, connection_list = create_dendrogram_map(radial, branches, layers)
#xs .+= rand(length(xs))/1000
#ys .+= rand(length(ys))/1000

connections = connection_matrix(connection_list, m = length(xs), n = length(ys))
dist_f(x) = 0.001 #constant distance function 
cell_map = CellMap(xs, ys, connections, distance_function = dist_f);
# Only do if there is GPU
#cell_map.strength = -cell_map.strength
cell_map_GPU = cell_map |> make_GPU
# [run the model]____________________________________________________________________________#
p0_dict = SAC_p0_dict()
p0_dict["g_ACh"] = 0.0
p0_dict["g_GABA"] = 0.0
p0_dict["g_GLUT"] = 1.0
p0_dict["I_app"] = 0.0
#p0_dict["g_K"] = LinRange(1.0, 10.0, size(cell_map)) |> collect
p0 = extract_p0(p0_dict)

u0_dict= SAC_u0_dict(n_cells = size(cell_map))
u0_dict["g"][1] = 2.0
u0 = extract_u0(u0_dict) 
u0 = u0 |> CuArray{Float32}
# 3) Define the problem
tspan = (0.0, 5e3)

f_PDE(du, u, p, t) = SAC_GAP(du, u, p, t, cell_map_GPU)
prob = SDEProblem(f_PDE, noise2D, u0, tspan, p0)

@time sol = solve(prob, 
    #SOSRI(), #This seems to be the best solver option
    SOSRA(),
    #dtmin = 0.001,
    reltol =  0.2, abstol = 2.0, 
    progress=true, progress_steps=1
);

# [Plot the model]___________________________________________________________#
CUDA.allowscalar(true) #allow GPU operations to be offloaded to CPU 
Time = LinRange(sol.t[1], sol.t[end], 5000)
It = hcat(map(t -> sol(t)[:,1], Time)...)|>Array
vt = hcat(map(t -> sol(t)[:,2], Time)...)|>Array
nt = hcat(map(t -> sol(t)[:,3], Time)...)|>Array
mt = hcat(map(t -> sol(t)[:,4], Time)...)|>Array
ht = hcat(map(t -> sol(t)[:,5], Time)...)|>Array
ct = hcat(map(t -> sol(t)[:,6], Time)...)|>Array
at = hcat(map(t -> sol(t)[:,7], Time)...)|>Array
bt = hcat(map(t -> sol(t)[:,8], Time)...)|>Array
et = hcat(map(t -> sol(t)[:,9], Time)...)|>Array
it = hcat(map(t -> sol(t)[:,10], Time)...)|>Array
Wt = hcat(map(t -> sol(t)[:,11], Time)...)|>Array

fig1 = Figure(size = (400,800))
ax1a = Axis(fig1[1,1], title = "Calcium Imageing", xlabel = "nx", ylabel = "ny")
ax1b = Axis(fig1[2,1], title = "Calcium ROIs", xlabel = "time (ms)", ylabel = "Ct")
ax1c = Axis(fig1[3,1], title = "Voltage Traces", xlabel = "time (ms)", ylabel = "Vt (mV)")
ax1d = Axis(fig1[4,1], title = "Current Traces", xlabel = "time (ms)", ylabel = "Vt (mV)")

rowsize!(fig1.layout, 1, Relative(1/2)) #Make the cell plot larger

rows, cols, data = findnz(cell_map.strength)
for i in eachindex(rows)
    r = rows[i]
    c = cols[i]
    d = data[i]
    lines!(ax1a, [xs[r], xs[c]], [ys[r], ys[c]], color = :black)#[d], colorrange = (minimum(data), maximum(data)), alpha = 0.1)
end
sctV = scatter!(ax1a, xs, ys, markersize = 15.0, color = sol(0.0)[:,2]|>Array, colorrange = (minimum(vt), maximum(vt)))
cell_map_GPU |> typeof

for n in 1:size(cell_map)
    lines!(ax1b, Time, ct[n, :])
    lines!(ax1c, Time, vt[n, :])
    lines!(ax1d, Time, It[n, :])
end
tickerb = vlines!(ax1b, [0.0])
tickerc = vlines!(ax1c, [0.0])

display(fig1)

n_frames = 1000
animate_t = LinRange(0.0, sol.t[end], n_frames)
dt = animate_t[2] - animate_t[1]
fps = round(Int64, (1/dt) * 1000)

GLMakie.record(fig1, "test/SAC_model_tests/data/branch_model_animation.mp4", animate_t, framerate = 8) do t
	println(t)
	c = sol(t)[:, 2] |> Array
	sctV.color = c
    tickerb[1] = [t]
    tickerc[1] = [t]
end   
using ElectroPhysiology, PhysiologyModeling
using PhysiologyPlotting, GLMakie

using CUDA #Running GPU in this section may be overkill
#CUDA.allowscalar(false)
using SparseArrays
using LinearAlgebra
using StatsBase, Statistics
using DiffEqCallbacks
import PhysiologyModeling.PresetTimeCallback
import PhysiologyModeling.create_dendrogram_map
import PhysiologyModeling.euclidean_distance
#%%=[Run branch generation]==================================================================#
radial = 5
branches = 2
layers = 5

xs1, ys1, connection_list1 = create_dendrogram_map(radial, branches, layers; origin = (0.0, 0.0))
xs2, ys2, connection_list_0 = create_dendrogram_map(radial, branches, layers; origin = (0.075,0.0))
x_max = size(xs1,1)
y_max = size(ys1,1)
connection_list2 = map(x -> (x[1]+x_max, x[2]+y_max, x[3]), connection_list_0)

xs = [xs1; xs2]
ys = [ys1; ys2]
connection_list = [connection_list1; connection_list2]
connections = connection_matrix(connection_list, m = length(xs), n = length(ys))
dist_f(p1, p2) = 10.0 #constant distance function 
V_MAP = CellMap(xs, ys, connections, distance_function = dist_f)# |> make_GPU;
#Check which connections are on the edge 
all_ys = map(i -> i[2], connection_list) 
all_ys_cm = countmap(all_ys)
all_ys_count = map(i -> all_ys_cm[i], all_ys)
no_release = all_ys[findall(all_ys_count .!= 1)]

#=[Run Acetylcholine propagation]==================================================================#
connection_list = connect_neighbors_radius(xs, ys, 0.01, self_connecting = false)
connections = connection_matrix(connection_list, m = length(xs), n = length(ys))
#Make the acetylcholine map
ACH_dist_func(p1, p2) = 0.05*euclidean_distance(p1, p2) 
E_MAP = CellMap(xs, ys, connections; distance_function = ACH_dist_func)# |> make_GPU;
#Make the GABA map between the two cells
GABA_dist_func(p1, p2) = 0.05*euclidean_distance(p1, p2)
I_MAP = CellMap(xs, ys, connections; distance_function = GABA_dist_func)# |> make_GPU;

#=[Simulate the stimulus]==================================================================#
p0_dict = SAC_p0_dict()
#p0_dict["g_ACh"] = 0.0
#p0_dict["g_GABA"] = 0.0
p0_dict["g_GLUT"] = 0.05
p0_dict["τq"] = 1500.0 
#Null out ACh release from all but edge dendrites
arr = fill(p0_dict["ρe"], size(xs,1))
map(i -> arr[i] = 0.0, no_release)
p0_dict["ρe"] = arr
#Null out GABA release from all but edge dendrites
arr = fill(p0_dict["ρi"], size(xs,1))
map(i -> arr[i] = 0.0, no_release)
p0_dict["ρi"] = arr

#Null out GABA release from all but edge dendrites
p0_dict["γg"] = fill(p0_dict["γg"], size(xs,1))
p0_dict["γg"][1] = 0.0
p0_dict["γg"][size(xs1, 1)+1] = 0.0

p0_dict["g_K"] = fill(4.0, size(xs,1))
p0_dict["g_K"][1] = 10.0
p0_dict["g_K"][size(xs1, 1)+1] = 10.0

p0 = extract_p0(p0_dict)

u0_dict= SAC_u0_dict(n_cells = size(xs,1))
u0 = extract_u0(u0_dict)# |> CuArray{Float32}

# 3) Define the problem
tspan = (0.0, 5000.0)
dt = 1.0
tseries = tspan[1]:dt:tspan[2]

f_PDE(du, u, p, t) = SAC_GAP(du, u, p, t, V_MAP, E_MAP, I_MAP)
prob = SDEProblem(f_PDE, noise2D, u0, tspan, p0)

#create a callback for the glutamate pulse
n_stops = 10
x_stops = LinRange(minimum(xs), maximum(xs), n_stops)

fn_affect!(integrator) = apply_glutamate_affect!(integrator; 
    xs = xs, 
    n_stops = n_stops, x_stops = x_stops,
    dt_pulse = 250.0, pulse_start = 500.0
)

cb = PresetTimeCallback(tseries, fn_affect!)
@time sol = solve(prob, 
    #SOSRI(), #This seems to be the best solver option
    SOSRA(),
    tstops = tseries, callback = cb,
    reltol =  0.2, abstol = 2.0, 
    progress=true, progress_steps=1
);

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
gt = hcat(map(t -> sol(t)[:,11], Time)...)|>Array
qt = hcat(map(t -> sol(t)[:,12], Time)...)|>Array
Wt = hcat(map(t -> sol(t)[:,13], Time)...)|>Array

#%%=[Plot the model]==================================================================#
save_fn = "Modeling/Results/branch_model_animation.mp4"

fig1 = Figure(size = (1000,700))
ax1a = Axis(fig1[1,1], title = "Calcium Imaging", xlabel = "nx", ylabel = "ny")
ax1b = Axis(fig1[2,1], title = "Calcium", xlabel = "time (ms)", ylabel = "Ct (mM)")
ax1c = Axis(fig1[3,1], title = "Voltage", xlabel = "time (ms)", ylabel = "vt (mV)")

ax1d = Axis(fig1[1,2], title = "GABA Imaging", xlabel = "nx", ylabel = "ny")
ax1e = Axis(fig1[2,2], title = "Current", xlabel = "time (ms)", ylabel = "It (pA)")
ax1f = Axis(fig1[3,2], title = "GABA", xlabel = "time (ms)", ylabel = "it (mM)")

ax1g = Axis(fig1[1,3], title = "Glutamate Imaging", xlabel = "nx", ylabel = "ny")
ax1h = Axis(fig1[2,3], title = "Glutamate", xlabel = "time (ms)", ylabel = "Gt (mM)")
ax1i = Axis(fig1[3,3], title = "G-protein", xlabel = "time (ms)", ylabel = "qt")

rowsize!(fig1.layout, 1, Relative(1/2)) #Make the cell plot larger
display(fig1)
rows, cols, data = findnz(V_MAP.strength)
for i in eachindex(rows)
    r = rows[i]
    c = cols[i]
    d = data[i]
    lines!(ax1a, [xs[r], xs[c]], [ys[r], ys[c]], color = :black)#[d], colorrange = (minimum(data), maximum(data)), alpha = 0.1)
    lines!(ax1d, [xs[r], xs[c]], [ys[r], ys[c]], color = :black)#[d], colorrange = (minimum(data), maximum(data)), alpha = 0.1)
    lines!(ax1g, [xs[r], xs[c]], [ys[r], ys[c]], color = :black)#[d], colorrange = (minimum(data), maximum(data)), alpha = 0.1)
end
sctC = scatter!(ax1a, xs, ys, markersize = 15.0, color = sol(0.0)[:,6]|>Array, colorrange = (minimum(ct), maximum(ct)))
sctI = scatter!(ax1d, xs, ys, markersize = 15.0, color = sol(0.0)[:,10]|>Array, colorrange = (minimum(it), maximum(it)))
sctG = scatter!(ax1g, xs, ys, markersize = 15.0, color = sol(0.0)[:,11]|>Array, colorrange = (minimum(gt), maximum(gt)))

for n in 1:size(xs,1)
    lines!(ax1b, Time, ct[n, :])
    lines!(ax1c, Time, vt[n, :])
    
    lines!(ax1e, Time, It[n, :])
    lines!(ax1f, Time, it[n, :])

    lines!(ax1h, Time, gt[n, :])
    lines!(ax1i, Time, qt[n, :])
end

lines!(ax1b, Time, ct[1, :], color = :red)
lines!(ax1c, Time, vt[1, :], color = :red)

lines!(ax1e, Time, It[1, :], color = :red)
lines!(ax1f, Time, it[1, :], color = :red)

lines!(ax1h, Time, gt[1, :], color = :red)
lines!(ax1i, Time, qt[1, :], color = :red)

tickerb = vlines!(ax1b, [0.0])
tickerc = vlines!(ax1c, [0.0])

tickere = vlines!(ax1e, [0.0])
tickerf = vlines!(ax1f, [0.0])

tickerh = vlines!(ax1h, [0.0])
tickeri = vlines!(ax1i, [0.0])

display(fig1)

#%%=[Animate]==================================================================#
n_frames = 500
animate_t = LinRange(0.0, sol.t[end], n_frames)
dt = animate_t[2] - animate_t[1]
fps = round(Int64, (1/dt) * 1000)

GLMakie.record(fig1, save_fn, animate_t, framerate = fps) do t
	println(t)
	sctC.color = sol(t)[:, 6] |> Array
    sctI.color = sol(t)[:, 10] |> Array
    sctG.color = sol(t)[:, 11] |> Array

    tickerb[1] = [t]
    tickerc[1] = [t]
    tickere[1] = [t]
    tickerf[1] = [t]
    tickerh[1] = [t]
    tickeri[1] = [t]
end   

#%%=[Plot the avergaes]==================================================================#
fig2 = Figure(size = (1000,700))
fig1.layout.width
ax1a = Axis(fig2[1,1], title = "Calcium Imaging", xlabel = "nx", ylabel = "ny")
ax1b = Axis(fig2[1,2], title = "GABA Imaging", xlabel = "nx", ylabel = "ny")
ax1c = Axis(fig2[1,3], title = "Glutamate Imaging", xlabel = "nx", ylabel = "ny")

rows, cols, data = findnz(V_MAP.strength)
for i in eachindex(rows)
    r = rows[i]
    c = cols[i]
    d = data[i]
    lines!(ax1a, [xs[r], xs[c]], [ys[r], ys[c]], color = :black)#[d], colorrange = (minimum(data), maximum(data)), alpha = 0.1)
    lines!(ax1b, [xs[r], xs[c]], [ys[r], ys[c]], color = :black)#[d], colorrange = (minimum(data), maximum(data)), alpha = 0.1)
    lines!(ax1c, [xs[r], xs[c]], [ys[r], ys[c]], color = :black)#[d], colorrange = (minimum(data), maximum(data)), alpha = 0.1)
end

c_avg = mean(ct, dims = 2)[:,1]
i_avg = mean(it, dims = 2)[:,1]
g_avg = mean(gt, dims = 2)[:,1]

sctC = scatter!(ax1a, xs, ys, markersize = 15.0, color = c_avg, colorrange = (minimum(c_avg), maximum(c_avg)))
sctI = scatter!(ax1b, xs, ys, markersize = 15.0, color = i_avg, colorrange = (minimum(i_avg), maximum(i_avg)))
sctG = scatter!(ax1c, xs, ys, markersize = 15.0, color = g_avg, colorrange = (minimum(g_avg), maximum(g_avg)))

display(fig2)
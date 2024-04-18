using Revise
using Pkg; Pkg.activate(".")
using ElectroPhysiology, PhysiologyModeling
Pkg.activate("test") #Activate the testing environment
using PhysiologyPlotting, GLMakie
#using CUDA #Running GPU in this section may be overkill
#CUDA.allowscalar(false)
using SparseArrays
using LinearAlgebra
using StatsBase, Statistics
save_fn = "test/SAC_model_tests/data/branch_model_animation.mp4"
import .PhysiologyModeling.DiffEqCallbacks
#%%=[Run branch generation]__________________________________________________________________________________#
radial = 5
branches = 2
layers = 5

xs, ys, connection_list = create_dendrogram_map(radial, branches, layers)
#xs .+= rand(length(xs))/1000
#ys .+= rand(length(ys))/1000

connections = connection_matrix(connection_list, m = length(xs), n = length(ys))
dist_f(x) = 10.0 #constant distance function 
cell_map = CellMap(xs, ys, connections, distance_function = dist_f);
# Only do if there is GPU
#cell_map.strength = -cell_map.strength
#cell_map = cell_map |> make_GPU



#%% [run the model]____________________________________________________________________________#
p0_dict = SAC_p0_dict()
p0_dict["g_ACh"] = 0.0
p0_dict["g_GABA"] = 0.0
p0_dict["g_GLUT"] = 1.0
p0_dict["g_K"] = [10.0, fill(4.0, size(cell_map)-1)...] 
p0 = extract_p0(p0_dict)

u0_dict= SAC_u0_dict(n_cells = size(cell_map))
u0 = extract_u0(u0_dict) 
#u0 = u0 |> CuArray{Float32}

# 3) Define the problem
tspan = (0.0, 5000.0)
dt = 1.0
tseries = tspan[1]:dt:tspan[2]

f_PDE(du, u, p, t) = SAC_GAP(du, u, p, t, cell_map)
prob = SDEProblem(f_PDE, noise2D, u0, tspan, p0)

#create a callback for the glutamate pulse
n_stops = 10
x_stops = LinRange(minimum(xs), maximum(xs), n_stops)

fn_affect!(integrator) = apply_glutamate_affect!(integrator; 
    xs = xs, 
    n_stops = n_stops, x_stops = x_stops,
    dt_pulse = 250.0, pulse_start = 500.0
)
import .PhysiologyModeling.PresetTimeCallback
cb = PresetTimeCallback(tseries, fn_affect!)
@time sol = solve(prob, 
    #SOSRI(), #This seems to be the best solver option
    SOSRA(),
    tstops = tseries, callback = cb,
    reltol =  0.2, abstol = 2.0, 
    progress=true, progress_steps=1
);

# 


# [Plot the model]___________________________________________________________#
#CUDA.allowscalar(true) #allow GPU operations to be offloaded to CPU 
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

#%% Figure 1
fig1 = Figure(size = (500,500))
ax1a = Axis(fig1[1,1], aspect = 1)

rows, cols, data = findnz(cell_map.strength)
for i in eachindex(rows)
    r = rows[i]
    c = cols[i]
    d = data[i]
    lines!(ax1a, [xs[r], xs[c]], [ys[r], ys[c]], color = :black)#[d], colorrange = (minimum(data), maximum(data)), alpha = 0.1)
end
texts = map(i -> "$i", axes(cell_map))
n_stops = 4
x_stops = LinRange(minimum(xs), maximum(xs), n_stops)

for i in 1:n_stops-1
    println(i)
    sctV = scatter!(ax1a, xs[findall(x_stops[i] .<= xs .<= x_stops[i+1])], ys[findall(x_stops[i] .<= xs .<= x_stops[i+1])], 
        markersize = 35.0, color = [i], colorrange = (1, n_stops), depth = 2 
    )
end
sctV = scatter!(ax1a, xs, ys, markersize = 25.0, color = :gold, depth = 2)

text!(ax1a, xs, ys, text = texts, align = (:center, :center))
display(fig1)

#%% Plot the map
fig2 = Figure(size = (1000,700))
fig2.layout.width
ax1a = Axis(fig2[1,1], title = "Calcium Imaging", xlabel = "nx", ylabel = "ny")
ax1b = Axis(fig2[2,1], title = "Calcium", xlabel = "time (ms)", ylabel = "Ct (mM)")
ax1c = Axis(fig2[3,1], title = "Voltage", xlabel = "time (ms)", ylabel = "vt (mV)")

ax1d = Axis(fig2[1,2], title = "GABA Imaging", xlabel = "nx", ylabel = "ny")
ax1e = Axis(fig2[2,2], title = "Current", xlabel = "time (ms)", ylabel = "It (pA)")
ax1f = Axis(fig2[3,2], title = "GABA", xlabel = "time (ms)", ylabel = "it (mM)")

ax1g = Axis(fig2[1,3], title = "Glutamate Imaging", xlabel = "nx", ylabel = "ny")
ax1h = Axis(fig2[2,3], title = "Glutamate", xlabel = "time (ms)", ylabel = "Gt (mM)")
ax1i = Axis(fig2[3,3], title = "G-protein", xlabel = "time (ms)", ylabel = "qt")

rowsize!(fig2.layout, 1, Relative(1/2)) #Make the cell plot larger
display(fig2)
rows, cols, data = findnz(cell_map.strength)
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

for n in 1:size(cell_map)
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

display(fig2)

n_frames = 500
animate_t = LinRange(0.0, sol.t[end], n_frames)
dt = animate_t[2] - animate_t[1]
fps = round(Int64, (1/dt) * 1000)

GLMakie.record(fig2, save_fn, animate_t, framerate = fps/10) do t
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
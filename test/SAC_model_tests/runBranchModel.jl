using Revise
using Pkg; Pkg.activate(".")
using ElectroPhysiology, PhysiologyModeling
Pkg.activate("test") #Activate the testing environment
using PhysiologyPlotting, GLMakie
#using CUDA #Running GPU in this section may be overkill
#CUDA.allowscalar(false)
using SparseArrays
using LinearAlgebra

save_fn = "test/SAC_model_tests/data/branch_model_animation.mp4"
import .PhysiologyModeling.DiffEqCallbacks
#%%=[Run branch generation]__________________________________________________________________________________#
radial = 4
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

fig = Figure(size = (500,500))
ax1a = Axis(fig[1,1], aspect = 1)

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
        markersize = 35.0, color = [i], colorrange = (1, n_stops), 
    )
end
sctV = scatter!(ax1a, xs, ys, markersize = 25.0, color = :gold)

text!(ax1a, xs, ys, text = texts, align = (:center, :center))
display(fig)

#%% [run the model]____________________________________________________________________________#
p0_dict = SAC_p0_dict()
p0_dict["g_ACh"] = 0.0
p0_dict["g_GABA"] = 0.0
p0_dict["g_GLUT"] = 1.0
#p0_dict["g_K"] = LinRange(1.0, 10.0, size(cell_map)) |> collect
p0 = extract_p0(p0_dict)

u0_dict= SAC_u0_dict(n_cells = size(cell_map))
#u0_dict["g"][5] = 1.0
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
    n_stops = n_stops, x_stops = x_stops,
    dt_pulse = 250.0, pulse_start = 500.0
)
cb = PresetTimeCallback(tseries, affect!)

@time sol = solve(prob, 
    #SOSRI(), #This seems to be the best solver option
    SOSRA(),
    tstops = tseries, callback = cb,
    reltol =  0.2, abstol = 2.0, 
    progress=true, progress_steps=1
);

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

fig1 = Figure(size = (1800,450))
ax1a = Axis(fig1[1:2,1], title = "Calcium Imageing", xlabel = "nx", ylabel = "ny")
ax1b = Axis(fig1[1,2], title = "Current", xlabel = "time (ms)", ylabel = "It (pA)")
ax1c = Axis(fig1[2,2], title = "Voltage", xlabel = "time (ms)", ylabel = "vt (mV)")

ax1d = Axis(fig1[1,3], title = "Calcium", xlabel = "time (ms)", ylabel = "Ct (mM)")
ax1e = Axis(fig1[2,3], title = "GABA", xlabel = "time (ms)", ylabel = "it (mM)")

ax1f = Axis(fig1[1,4], title = "Glutamate", xlabel = "time (ms)", ylabel = "Gt (mM)")
ax1g = Axis(fig1[2,4], title = "G-protein", xlabel = "time (ms)", ylabel = "qt")

colsize!(fig1.layout, 1, Relative(1/4)) #Make the cell plot larger
rows, cols, data = findnz(cell_map.strength)
for i in eachindex(rows)
    r = rows[i]
    c = cols[i]
    d = data[i]
    lines!(ax1a, [xs[r], xs[c]], [ys[r], ys[c]], color = :black)#[d], colorrange = (minimum(data), maximum(data)), alpha = 0.1)
end
sctV = scatter!(ax1a, xs, ys, markersize = 15.0, color = sol(0.0)[:,6]|>Array, colorrange = (minimum(ct), maximum(ct)))

for n in 1:size(cell_map)
    lines!(ax1b, Time, It[n, :])
    lines!(ax1c, Time, vt[n, :])
    lines!(ax1d, Time, ct[n, :])
    lines!(ax1e, Time, it[n, :])
    lines!(ax1f, Time, gt[n, :])
    lines!(ax1g, Time, qt[n, :])
end
tickerb = vlines!(ax1b, [0.0])
tickerc = vlines!(ax1c, [0.0])
tickerd = vlines!(ax1d, [0.0])
tickere = vlines!(ax1e, [0.0])
tickerf = vlines!(ax1f, [0.0])
tickerg = vlines!(ax1g, [0.0])

display(fig1)


n_frames = 500
animate_t = LinRange(0.0, sol.t[end], n_frames)
dt = animate_t[2] - animate_t[1]
fps = round(Int64, (1/dt) * 1000)

GLMakie.record(fig1, save_fn, animate_t, framerate = fps/10) do t
	println(t)
	c = sol(t)[:, 6] |> Array
	sctV.color = c
    tickerb[1] = [t]
    tickerc[1] = [t]
    tickerd[1] = [t]
    tickere[1] = [t]
    tickerf[1] = [t]
    tickerg[1] = [t]
end   
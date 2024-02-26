using Revise 
using ElectroPhysiology
using PhysiologyModeling
import PhysiologyModeling: CVODE_BDF, ring
using Pkg; Pkg.activate("test") #Activate the testing environment
using PhysiologyPlotting, GLMakie
using CUDA
CUDA.allowscalar(false)
using SparseArrays
using LinearAlgebra

#%%

#%% 1) determine the domains and spacing of cells. 
domain_x = (xmin, xmax) = (0.0, 1.0)
domain_y = (ymin, ymax) = (0.0, 1.0)
dx = dy = 0.05 #Mean distribution is 40-50 micron (WR taylor et al)

#2) create a random distribution of cells and their radii
ncells = 200
cells = rand(xmin:dx:xmax, ncells, 2)
xs = cells[:, 1]
ys = cells[:, 2]
#%%
#cells = even_map(xmin = xmin, dx = dx, xmax = xmax, ymin = ymin, dy = dy, ymax = ymax)
radii = fill(0.18, size(cells, 1)) #Switch this on to get constant radii
dist_func(d) = ring(d; max_strength = 0.1, max_dist = 0.15, slope = 0.05)
cell_map = CellMap(cells, radii; distance_function = dist_func);
#cell_map = cell_map |> make_GPU

p0_dict = SAC_p0_dict()
p0_dict["g_ACh"] = 0.215
p0_dict["g_GABA"] = 0.0
p0_dict["g_W"] = 0.0 
p0_dict["I_app"] = -20.0
p0 = extract_p0(p0_dict) 

u0_dict= SAC_u0_dict(mode = :PDE, ncells = ncells)
u0_dict["v"][10] = 20.0
u0 = extract_u0(u0_dict)# |> CuArray{Float32}

#3) Define the problem
tspan = (0.0, 100.0)
f_PDE(du, u, p, t) = SAC_PDE(du, u, p, t, cell_map)
prob = SDEProblem(f_PDE, noise2D, u0, tspan, p0)
@time sol = solve(prob, 
     SOSRI(), #This seems to be the best solver option
      force_dtmin = true, 
     progress=true, progress_steps=1)
#save("data.jld", "initial_cond", sol[end])

#%% Start plotting
 CUDA.allowscalar(true) #allow GPU operations to be offloaded to CPU 
Time = sol.t[1]:10:sol.t[end]
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

dist_func(0.1)
#Create the plot
fig1 = Figure(size = (400,800))
ax1a = Axis(fig1[1,1], title = "Voltage")
ax1b = Axis(fig1[2,1], title = "Voltage Traces")
sctV = scatter!(ax1a, cells[:,1], cells[:,2], color = sol(0.0)[:,1]|>Array, colorrange = (minimum(vt), maximum(vt)), markersize = 20.0)
for n in 1:ncells
     lines!(ax1b, Time, vt[n, :])
end
ticker = vlines!(ax1b, [0.0])
for (i,c) in enumerate(connected_idxs[1])
     xs_origin = cell_map.xs[idx]
     ys_origin = cell_map.ys[idx]
     xs_end = cell_map.xs[c]
     ys_end = cell_map.ys[c]

     scatter!(ax1a, xs_origin, ys_origin, markersize = 0.4, alpha = 0.1, markerspace = :data)
     scatter!(ax1a, xs_origin, ys_origin, color = :red)
     lines!(ax1a, [xs_origin, xs_end], [ys_origin, ys_end], colormap = :viridis, color = cell_map.strength[idx, c], colorrange = (0.0, 0.2))
end
display(fig1)

n_frames = 1000
animate_t = LinRange(0.0, sol.t[end], n_frames)
dt = animate_t[2] - animate_t[1]
fps = round(Int64, (1/dt) * 1000)

GLMakie.record(fig1, "test/SAC_model_tests/model_animation.mp4", animate_t, framerate = 8) do t
	println(t)
	v = sol(t)[:, 2] |> Array
	sctV.color = v
     ticker[1] = [t]
end


#%%
ax2 = Axis(fDIFF[1,2]; aspect=(1, 1, 1), title = "Acetylcholine")
ax3 = Axis(fDIFF[1,3]; aspect=(1, 1, 1), title = "GABA")

ax4 = Axis(fDIFF[2,1]; aspect=(1, 1, 1), title = "Calcium")
ax5 = Axis(fDIFF[2,2]; aspect=(1, 1, 1), title = "cAMP")
ax6 = Axis(fDIFF[2,3]; aspect=(1, 1, 1), title = "TREK1")

surf1 = surface!(ax1, cells[:, 1], cells[:, 2], sol(0.0)[:, 1])
surf2 = surface!(ax2, cells[:, 1], cells[:, 2], sol(0.0)[:, 8])
surf3 = surface!(ax3, cells[:, 1], cells[:, 2], sol(0.0)[:, 9])

surf4 = surface!(ax4, cells[:, 1], cells[:, 2], sol(0.0)[:, 5])
surf5 = surface!(ax5, cells[:, 1], cells[:, 2], sol(0.0)[:, 6])
surf6 = surface!(ax6, cells[:, 1], cells[:, 2], sol(0.0)[:, 7])

zlims!(ax1, zlims1)
zlims!(ax2, zlims2)
zlims!(ax3, zlims3)
zlims!(ax4, zlims4)
zlims!(ax5, zlims5)
zlims!(ax6, zlims6)

#display(fDIFF)
n_frames = 1000
animate_t = LinRange(0.0, sol.t[end], n_frames)
dt = animate_t[2] - animate_t[1]
fps = round(Int64, (1/dt) * 1000)

GLMakie.record(fDIFF, "test/SAC_model_tests/model_animation.mp4", animate_t, framerate = 8) do t
	println(t)
	v = sol(t)[:, 1]
     e = sol(t)[:, 8]
     i = sol(t)[:, 9]
     c = sol(t)[:, 5]
     a = sol(t)[:, 6]
     b = sol(t)[:, 7]
	surf1[3] = v
     surf2[3] = e
     surf3[3] = i  
     surf4[3] = c
     surf5[3] = a
     surf6[3] = b
end

#%%
CUDA.allowscalar(true) #
fSDE = Figure(size = (1800, 800))
ax1 = Axis(fSDE[1,1], title = "Voltage (Vt)")
ax2 = Axis(fSDE[2,1], title = "K Repol. (Nt)")
ax3 = Axis(fSDE[3,1], title = "Na Gating (Mt)")
ax4 = Axis(fSDE[4,1], title = "Na Close (Ht)")

ax5 = Axis(fSDE[1,2], title = "Calcium (Ct)")
ax6 = Axis(fSDE[2,2], title = "cAMP (At)")
ax7 = Axis(fSDE[3,2], title = "TREK (Bt)")

ax8 = Axis(fSDE[1,3], title = "ACh (Et)")
ax9 = Axis(fSDE[2,3], title = "GABA (It)")

ax10 = Axis(fSDE[1,4], title = "Noise (Wt)")
Time = LinRange(sol.t[1], sol.t[end], 1000)
for i in rand(1:size(sol, 1), 10)
     lines!(ax1, Time, map(t -> sol(t)[i, 1], Time))
     lines!(ax2, Time, map(t -> sol(t)[i, 2], Time))
     lines!(ax3, Time, map(t -> sol(t)[i, 3], Time))
     lines!(ax4, Time, map(t -> sol(t)[i, 4], Time))
     lines!(ax5, Time, map(t -> sol(t)[i, 5], Time))
     lines!(ax6, Time, map(t -> sol(t)[i, 6], Time))
     lines!(ax7, Time, map(t -> sol(t)[i, 7], Time))
     lines!(ax8, Time, map(t -> sol(t)[i, 8], Time))
     lines!(ax9, Time, map(t -> sol(t)[i, 9], Time))
     lines!(ax10, Time, map(t -> sol(t)[i, 10], Time))
end
display(fSDE)
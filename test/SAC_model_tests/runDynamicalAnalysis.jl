using Revise
using ElectroPhysiology, PhysiologyModeling
using Pkg; Pkg.activate("test") #Activate the testing environment
using PhysiologyPlotting, GLMakie
using BenchmarkTools

#%% Create the figure basis
fig1 = Figure(size = (1600, 800))
ga = fig1[1, 1] = GridLayout()
axA1 = Axis(ga[1, 1], xlabel = "Time (ms)", ylabel = "Voltage (Vt)")
axA2 = Axis(ga[2, 1], xlabel = "Time (ms)", ylabel = "Repol (Nt)")

gb = fig[1, 2] = GridLayout()
axB1 = Axis3(gb[1,1], title = "K Repolarization (Nt)", xlabel = "I_app", ylabel = "Voltage. (Vt) ", zlabel = "K Repol. (N)")

# Set some initial parameters
u0_dict = SAC_u0_dict()
p0_dict = SAC_p0_dict(keyset = :dynamical_analysis)
u0_dict["v"] = -65.0

# After setting the inital parameters, extract them
u0 = extract_u0(u0_dict)
p0 = extract_p0(p0_dict)
prob = ODEProblem(SAC_ODE, u0, (0.0, 300.0), p0) #Create the problem

#Set up the ensemble solution
n_traces = 100 #Number of traces
p_rng = LinRange(-20, 250, n_traces) #Set the parameter range
function prob_func(prob, i, repeat; idx = 1)
     pI = prob.p
     pI[idx] = p_rng[i]
     remake(prob, p = pI)
end
ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)

#Solve the ensemble solution
@time sim = solve(ensemble_prob, EnsembleDistributed(), trajectories = n_traces, progress = true, progress_steps = 1, maxiters = 1e7);

# Plot the results
for (i, sol) in enumerate(sim)
     println(i)
     x = p_rng[i]
     dt = 0.01
     Time = sol.t[1]:dt:sol.t[end]
     vt = map(t -> sol(t)[2], Time)
     nt = map(t -> sol(t)[3], Time)
     ct = map(t -> sol(t)[6], Time)

     lines!(axA1, Time, vt, colormap = :rainbow, color = [x], colorrange = (p_rng[1], p_rng[end]))
     lines!(axA2, Time, nt, colormap = :rainbow, color = [x], colorrange = (p_rng[1], p_rng[end]))
     lines!(axB1, [x], vt, nt, colormap = :rainbow, color = [x], colorrange = (p_rng[1], p_rng[end]))
end
display(fig1)

#%% Phase plane analysis
fig2 = Figure(size = (1200, 400))
ax1 = Axis3(fig2[1, 1], title = "Vt vector", xlabel = "Voltage (mV)", ylabel = "Repol (Vt)")
ax2 = Axis3(fig2[1, 2], title = "Nt vector", xlabel = "Voltage (mV)", ylabel = "Repol (Nt)")
ax3 = Axis(fig2[1, 3], title = "Nt vector", xlabel = "Voltage (mV)", ylabel = "Repol (Nt)")

p0_dict = SAC_p0_dict(keyset = :dynamical_analysis)
p0_dict["I_app"] = 10.0
u0 = extract_u0(u0_dict)
p0 = extract_p0(p0_dict)
prob = ODEProblem(SAC_ODE, u0, (0.0, 300.0), p0) #Create the problem

nx = ny = 50
xmap = LinRange(-70.0, 10.0, nx)
ymap = LinRange(-0.1, 1.1, ny)
phase_map = phase_plane(prob, xmap, ymap)
vt_vector = phase_map[:,:,1]
nt_vector = phase_map[:,:,2]
strength = sqrt.(vt_vector .^ 2 .+ nt_vector .^ 2) 

surface!(ax1, xmap, ymap, vt_vector)
surface!(ax2, xmap, ymap, nt_vector)

heatmap!(ax3, xmap, ymap, strength, alpha = 0.2)
arrows!(ax3, xmap, ymap, vt_vector, nt_vector, 
     arrowsize = 10.0, 
     lengthscale = 0.01, normalize = true,
     arrowcolor = strength|>vec, 
)

# finding nullclines
vt_nc, nt_nc = find_nullclines(prob, xmap, ymap)
vt_nc
nt_nc
lines!(ax1, xmap, nt_nc, color = :black, linewidth = 5.0)
lines!(ax2, vt_nc, ymap, color = :black, linewidth = 5.0)
display(fig2)

#%% Can we come up with a clever way to find the nullcline
#find it from the xperspective
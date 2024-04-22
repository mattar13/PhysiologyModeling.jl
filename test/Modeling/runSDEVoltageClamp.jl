#=[Import packages]============================================================#
using ElectroPhysiology, PhysiologyModeling
using PhysiologyPlotting, GLMakie
using DiffEqCallbacks

#%%=[Solving a single SDE for tspan]==========================================#
save_fn = "test/Modeling/Results/VC hold -65.png"
#Extract and modify initial conditions
u0_dict = SAC_u0_dict()
SAC_u0_dict
u0 = extract_u0(u0_dict)

#Specify the timespan
tspan = (0.0, 3e3)

#Set the stimulus parameters
p0_dict = SAC_p0_dict() #Extract parameters
p0_dict["g_GABA"] = 0.0 #Null GABA conductance
p0_dict["g_ACh"] = 0.0 #Null ACh conductance

p0_dict["stim_start"] = 200.0
p0_dict["VC"] = -15.0
p0_dict["stim_stop"] = 500.0

p0 = extract_p0(p0_dict)

#Set up the problem
fSDE_VC(du, u, p, t) = SAC_ODE_VC(du, u, p, t; k = 2000.0)
prob = SDEProblem(fSDE_VC, noise1D, u0, tspan, p0)
@time sol = solve(prob, SOSRI(), reltol = 2e-2, abstol = 2e-2, progress = true, progress_steps = 1)

# [Plot the solution]_________________________________________________________________________________________________________#
fig1 = Figure(size = (1000, 800))
ax1a = Axis(fig1[1,1], title = "I_ext (pA)")
ax2a = Axis(fig1[2,1], title = "Voltage (Vt)")
ax3a = Axis(fig1[3,1], title = "Noise (Wt)")
#=
ax1b = Axis(fig1[1,2], title = "K Repol. (Nt)")
ax2b = Axis(fig1[2,2], title = "Na Gating (Mt)")
ax3b = Axis(fig1[3,2], title = "Na Close (Ht)")

ax1c = Axis(fig1[1,3], title = "Calcium (Ct)")
ax2c = Axis(fig1[2,3], title = "cAMP (At)")
ax3c = Axis(fig1[3,3], title = "TREK (Bt)")

ax1d = Axis(fig1[1,4], title = "ACh (Et)")
ax2d = Axis(fig1[2,4], title = "GABA (It)")
ax3d = Axis(fig1[3,4], title = "Glutamate (Gt)")
=#
Time = sol.t
lines!(ax1a, Time, map(t -> sol(t)[1], Time))
lines!(ax2a, Time, map(t -> sol(t)[2], Time))
lines!(ax3a, Time, map(t -> sol(t)[13], Time))
#=
lines!(ax1b, Time, map(t -> sol(t)[3], Time))
lines!(ax2b, Time, map(t -> sol(t)[4], Time))
lines!(ax3b, Time, map(t -> sol(t)[5], Time))

lines!(ax1c, Time, map(t -> sol(t)[6], Time))
lines!(ax2c, Time, map(t -> sol(t)[7], Time))
lines!(ax3c, Time, map(t -> sol(t)[8], Time))

lines!(ax1d, Time, map(t -> sol(t)[9], Time))
lines!(ax2d, Time, map(t -> sol(t)[10], Time))
lines!(ax3d, Time, map(t -> sol(t)[11], Time))
=#
display(fig1)
save(save_fn, fig1)

#%% [Simulate a episodic current clamp experiment]_______________________________________________________________________________________#

#Specify the timespan
tspan = (0.0, 3e3)
#Extract and modify parameters
p0_dict = SAC_p0_dict()
p0_dict["g_GABA"] = 0.0
p0_dict["g_ACh"] = 0.0

p0_dict["stim_start"] = 200.0
p0_dict["VC"] = -10.0
p0_dict["stim_stop"] = 500.0

p0 = extract_p0(p0_dict)

#Extract and modify initial conditions
u0_dict = SAC_u0_dict()
u0 = extract_u0(u0_dict)

#Set up the problem
idx = par_idx("VC")
n_traces = 10
param_rng = LinRange(-100.0, 100, n_traces)

fSAC(du, u, p, t) = SAC_ODE_VC(du, u, p, t; hold = -95.0)
prob = SDEProblem(fSAC, noise1D, u0, tspan, p0)
function fSAC_ensemble(prob, i, repeat; idx = idx)
     pI = prob.p
     pI[idx] = param_rng[i]
     remake(prob, p = pI) 
end

ensemble_prob = EnsembleProblem(prob, prob_func = fSAC_ensemble)
@time sim = solve(ensemble_prob, SOSRI(), EnsembleDistributed(), trajectories = n_traces, 
     progress = true, progress_steps = 1, reltol = 0.01, abstol = 0.01, maxiters = 1e7)

# [Plot the solution]_____________________________________________________________________________________#
fig2 = Figure(size = (1000, 800))
ax1a = Axis(fig2[1,1], title = "I_ext (pA)")
ax2a = Axis(fig2[2,1], title = "Voltage (Vt)")
ax3a = Axis(fig2[3,1], title = "Noise (Wt)")

#ax1b = Axis(fig2[1,2], title = "K Repol. (Nt)")
#ax2b = Axis(fig2[2,2], title = "Na Gating (Mt)")
#ax3b = Axis(fig2[3,2], title = "Na Close (Ht)")

#ax1c = Axis(fig2[1,3], title = "Calcium (Ct)")
#ax2c = Axis(fig2[2,3], title = "cAMP (At)")
#ax3c = Axis(fig2[3,3], title = "TREK (Bt)")

#ax1d = Axis(fig2[1,4], title = "ACh (Et)")
#ax2d = Axis(fig2[2,4], title = "GABA (It)")
#ax3d = Axis(fig2[3,4], title = "Glutamate (Gt)")

Colorbar(fIC[5,1], limits = (param_rng[1], param_rng[end]), colormap = :viridis,
    vertical = false)
for (i, sol) in enumerate(sim)
     println(sol.t[end])
     Time = LinRange(sol.t[1], sol.t[end], 2000)
     x = param_rng[i]
     Time = sol.t
     Time = sol.t
     lines!(ax1a, Time, map(t -> sol(t)[1], Time))
     lines!(ax2a, Time, map(t -> sol(t)[2], Time))
     lines!(ax3a, Time, map(t -> sol(t)[13], Time))
     
     #lines!(ax1b, Time, map(t -> sol(t)[3], Time))
     #lines!(ax2b, Time, map(t -> sol(t)[4], Time))
     #lines!(ax3b, Time, map(t -> sol(t)[5], Time))
     
     #lines!(ax1c, Time, map(t -> sol(t)[6], Time))
     #lines!(ax2c, Time, map(t -> sol(t)[7], Time))
     #lines!(ax3c, Time, map(t -> sol(t)[8], Time))
     
     #lines!(ax1d, Time, map(t -> sol(t)[9], Time))
     #lines!(ax2d, Time, map(t -> sol(t)[10], Time))
     #lines!(ax3d, Time, map(t -> sol(t)[11], Time))
end

display(fig2)
using Pkg; Pkg.activate(".")
using Revise
using ElectroPhysiology
using PhysiologyModeling
using DiffEqCallbacks
Pkg.activate("test")
using PhysiologyPlotting
using GLMakie

#%% [Solving a single SDE for tspan]_______________________________________________________________________________________________#
tspan = (0.0, 300e3)
#Extract and modify parameters
p0_dict = SAC_p0_dict()
p0_dict["g_GABA"] = 0.0
p0_dict["g_ACh"] = 0.0
#Set the stimulus parameters
p0 = extract_p0(p0_dict)
#Extract and modify initial conditions
u0_dict = SAC_u0_dict()
u0 = extract_u0(u0_dict)
#Set up the problem
prob = SDEProblem(SAC_ODE, noise1D, u0, tspan, p0)
@time sol = solve(prob, SOSRI(), reltol = 2e-2, abstol = 2e-2, progress = true, progress_steps = 1)

# [Plot the solution]_________________________________________________________________________________________________________#
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
ax10 = Axis(fSDE[3,3], title = "Glutamate (Gt)")
ax11 = Axis(fSDE[4,3], title = "Gq (Qt)")

ax12 = Axis(fSDE[1,4], title = "Noise (Wt)")
ax13 = Axis(fSDE[2,4], title = "I_ext (pA)")

Time = sol.t
lines!(ax1, Time, map(t -> sol(t)[2], Time))
lines!(ax2, Time, map(t -> sol(t)[3], Time))
lines!(ax3, Time, map(t -> sol(t)[4], Time))
lines!(ax4, Time, map(t -> sol(t)[5], Time))
lines!(ax5, Time, map(t -> sol(t)[6], Time))
lines!(ax6, Time, map(t -> sol(t)[7], Time))
lines!(ax7, Time, map(t -> sol(t)[8], Time))
lines!(ax8, Time, map(t -> sol(t)[9], Time))
lines!(ax9, Time, map(t -> sol(t)[10], Time))
lines!(ax10, Time, map(t -> sol(t)[11], Time))
lines!(ax11, Time, map(t -> sol(t)[12], Time))
lines!(ax12, Time, map(t -> sol(t)[13], Time))
lines!(ax13, Time, map(t -> sol(t)[1], Time))

display(fSDE)
save("test/SAC_model_tests/data/SDESol_bKV.png", fSDE)
#%% [Simulate a current clamp experiment]_______________________________________________________________________________________#

#Specify the timespan
tspan = (0.0, 3e3)
#Extract and modify parameters
p0_dict = SAC_p0_dict()
p0_dict["g_GABA"] = 0.0
p0_dict["g_ACh"] = 0.0

p0_dict["stim_start"] = 200.0
p0_dict["I_app"] = -10.0
p0_dict["stim_stop"] = 500.0

p0 = extract_p0(p0_dict)

#Extract and modify initial conditions
u0_dict = SAC_u0_dict()
u0 = extract_u0(u0_dict)

#Set up the problem
idx = par_idx("I_app")
n_traces = 10
param_rng = LinRange(-20.0, 200, n_traces)

fSAC(du, u, p, t) = SAC_ODE_IC(du, u, p, t)
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
fIC = Figure(size = (1800, 800))
ax1 = Axis(fIC[1,1], title = "Voltage (Vt)")
ax2 = Axis(fIC[2,1], title = "K Repol. (Nt)")
ax3 = Axis(fIC[3,1], title = "Na Gating (Mt)")
ax4 = Axis(fIC[4,1], title = "Na Close (Ht)")

ax5 = Axis(fIC[1,2], title = "Calcium (Ct)")
ax6 = Axis(fIC[2,2], title = "cAMP (At)")
ax7 = Axis(fIC[3,2], title = "TREK (Bt)")

ax8 = Axis(fIC[1,3], title = "ACh (Et)")
ax9 = Axis(fIC[2,3], title = "GABA (It)")
ax10 = Axis(fIC[3,3], title = "Glutamate (Gt)")
ax11 = Axis(fIC[4,3], title = "Gq (Qt)")

ax12 = Axis(fIC[1,4], title = "Noise (Wt)")
ax13 = Axis(fIC[2,4], title = "I_ext (pA)")

Colorbar(fIC[5,1], limits = (param_rng[1], param_rng[end]), colormap = :viridis,
    vertical = false)
for (i, sol) in enumerate(sim)
     println(sol.t[end])
     Time = LinRange(sol.t[1], sol.t[end], 2000)
     x = param_rng[i]
     Time = sol.t
     ln1 = lines!(ax1, Time, map(t -> sol(t)[2], Time), color = x, colormap = :viridis, colorrange = (param_rng[1], param_rng[end]))
     lines!(ax2, Time, map(t -> sol(t)[3], Time), color = x, colormap = :viridis, colorrange = (param_rng[1], param_rng[end]))
     lines!(ax3, Time, map(t -> sol(t)[4], Time), color = x, colormap = :viridis, colorrange = (param_rng[1], param_rng[end]))
     lines!(ax4, Time, map(t -> sol(t)[5], Time), color = x, colormap = :viridis, colorrange = (param_rng[1], param_rng[end]))
     lines!(ax5, Time, map(t -> sol(t)[6], Time), color = x, colormap = :viridis, colorrange = (param_rng[1], param_rng[end]))
     lines!(ax6, Time, map(t -> sol(t)[7], Time), color = x, colormap = :viridis, colorrange = (param_rng[1], param_rng[end]))
     lines!(ax7, Time, map(t -> sol(t)[8], Time), color = x, colormap = :viridis, colorrange = (param_rng[1], param_rng[end]))
     lines!(ax8, Time, map(t -> sol(t)[9], Time), color = x, colormap = :viridis, colorrange = (param_rng[1], param_rng[end]))
     lines!(ax9, Time, map(t -> sol(t)[10], Time), color = x, colormap = :viridis, colorrange = (param_rng[1], param_rng[end]))
     lines!(ax10, Time, map(t -> sol(t)[11], Time), color = x, colormap = :viridis, colorrange = (param_rng[1], param_rng[end]))
     lines!(ax11, Time, map(t -> sol(t)[12], Time), color = x, colormap = :viridis, colorrange = (param_rng[1], param_rng[end]))
     lines!(ax12, Time, map(t -> sol(t)[13], Time), color = x, colormap = :viridis, colorrange = (param_rng[1], param_rng[end]))
     lines!(ax13, Time, map(t -> sol(t)[1], Time), color = x, colormap = :viridis, colorrange = (param_rng[1], param_rng[end]))
end

display(fIC)
 
#%% [Solving a single SDE by adding mGluR6 into the equation]______________________________________________________________________#
tspan = (0.0, 10e3)
dt = 1.0
tseries = tspan[1]:dt:tspan[2]

#Extract and modify parameters
p0_dict = SAC_p0_dict()
p0_dict["VC"] = -65.0
p0_dict["g_GABA"] = 0.0
p0_dict["g_ACh"] = 0.0
p0 = extract_p0(p0_dict)
#Extract and modify initial conditions
u0_dict = SAC_u0_dict()
u0_dict["g"] = 0.0
u0 = extract_u0(u0_dict)

#create a callback for the glutamate pulse
function affect!(integrator; t1 = 0.0, dtpulse = 2500.0, tend = 10e3)
     gt = 0.0
     for t0 in t1:dtpulse:tend
          gt += gauss_pulse(integrator.t; t0 = t0, peak_amp = 0.05)
     end
     integrator.u[11] = gt
end
cb = PresetTimeCallback(tseries, affect!)

#Set up the problem
prob_func(du, u, p, t) = SAC_ODE_VC(du, u, p, t)
prob = SDEProblem(prob_func, noise1D, u0, tspan, p0)
@time sol = solve(prob, SOSRI(), tstops = tseries, callback = cb, reltol = 2e-2, abstol = 2e-2, progress = true, progress_steps = 1)

# [Plot the solution]_______________________________________________________________________________________________________________#
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
ax10 = Axis(fSDE[3,3], title = "Glutamate (Gt)")
ax11 = Axis(fSDE[4,3], title = "Gq (Qt)")

ax12 = Axis(fSDE[1,4], title = "Noise (Wt)")
ax13 = Axis(fSDE[2,4], title = "I_ext (pA)")

Time = sol.t
lines!(ax1, Time, map(t -> sol(t)[2], Time))
lines!(ax2, Time, map(t -> sol(t)[3], Time))
lines!(ax3, Time, map(t -> sol(t)[4], Time))
lines!(ax4, Time, map(t -> sol(t)[5], Time))
lines!(ax5, Time, map(t -> sol(t)[6], Time))
lines!(ax6, Time, map(t -> sol(t)[7], Time))
lines!(ax7, Time, map(t -> sol(t)[8], Time))
lines!(ax8, Time, map(t -> sol(t)[9], Time))
lines!(ax9, Time, map(t -> sol(t)[10], Time))
lines!(ax10, Time, map(t -> sol(t)[11], Time))
lines!(ax11, Time, map(t -> sol(t)[12], Time))
lines!(ax12, Time, map(t -> sol(t)[13], Time))
lines!(ax13, Time, map(t -> sol(t)[1], Time))

display(fSDE)
save("test/SAC_model_tests/data/SDESol_GLUT.png", fSDE)
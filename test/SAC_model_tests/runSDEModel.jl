using PhysiologyModeling

using Pkg; Pkg.activate("test")
using PhysiologyPlotting
using GLMakie

#%% Run the example after here
tspan = (0.0, 4000.0)
SAC_p0_dict["I_app"] = 5.0
SAC_p0_dict["g_GABA"] = 0.0
SAC_p0_dict["g_ACh"] = 0.0
SAC_p0_dict["g_W"] = 0.0
SAC_p0_dict["g_TREK"] = 3.0
p0 = extract_p0(SAC_p0_dict)
u0 = extract_u0(SAC_u0_dict)
prob = SDEProblem(SAC_ODE_STIM, noise1D, u0, tspan, p0)
@time sol = solve(prob, SOSRI(), reltol = 0.01, abstol = 0.01, progress = true, progress_steps = 1)

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
Time = sol.t
lines!(ax1, Time, map(t -> sol(t)[1], Time))
lines!(ax2, Time, map(t -> sol(t)[2], Time))
lines!(ax3, Time, map(t -> sol(t)[3], Time))
lines!(ax4, Time, map(t -> sol(t)[4], Time))
lines!(ax5, Time, map(t -> sol(t)[5], Time))
lines!(ax6, Time, map(t -> sol(t)[6], Time))
lines!(ax7, Time, map(t -> sol(t)[7], Time))
lines!(ax8, Time, map(t -> sol(t)[8], Time))
lines!(ax9, Time, map(t -> sol(t)[9], Time))
lines!(ax10, Time, map(t -> sol(t)[10], Time))
display(fSDE)
#save("test/SAC_model_tests/SDESol.png", fSDE)


#%% Run an ensemble solution
tspan = (0.0, 300e3)
SAC_p0_dict["I_app"] = 5.0
SAC_p0_dict["g_GABA"] = 0.0
SAC_p0_dict["g_ACh"] = 0.0
SAC_p0_dict["g_W"] = 0.1
SAC_p0_dict["g_TREK"] = 4.0
p0 = extract_p0(SAC_p0_dict)
u0 = extract_u0(SAC_u0_dict)
prob = SDEProblem(SAC_ODE_STIM, noise1D, u0, tspan, p0)

n_traces = 20
initial_conditions = LinRange(-250, 100.0, n_traces)
function prob_func(prob, i, repeat; idx = 1)
     pI = prob.p
     pI[idx] = initial_conditions[i]
     println(pI)
     remake(prob, p = pI)
end

ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
@time sim = solve(ensemble_prob, SOSRI(), EnsembleDistributed(), trajectories = n_traces, 
     progress = true, progress_steps = 1, reltol = 0.01, abstol = 0.01, maxiters = 1e7)
#%%
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
for sol in sim
     println(sol.t[end])
     Time = sol.t[1]:1.0:4000.0
     lines!(ax1, Time, map(t -> sol(t)[1], Time))
     lines!(ax2, Time, map(t -> sol(t)[2], Time))
     lines!(ax3, Time, map(t -> sol(t)[3], Time))
     lines!(ax4, Time, map(t -> sol(t)[4], Time))
     lines!(ax5, Time, map(t -> sol(t)[5], Time))
     lines!(ax6, Time, map(t -> sol(t)[6], Time))
     lines!(ax7, Time, map(t -> sol(t)[7], Time))
     lines!(ax8, Time, map(t -> sol(t)[8], Time))
     lines!(ax9, Time, map(t -> sol(t)[9], Time))
     lines!(ax10, Time, map(t -> sol(t)[10], Time))
end
display(fSDE)
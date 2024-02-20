using PhysiologyModeling
using Pkg; Pkg.activate("test")
using PhysiologyPlotting
using GLMakie
#%% Run a Current clamp experiment
tspan = (0.0, 1e3)
dt = 1.0
tstops = tspan[1]:dt:tspan[end]
gE = rand(-200.0:0.0, length(tstops))
gI = rand(0.0:200.0, length(tstops))
p0_dict = SAC_p0_dict()
p0 = extract_p0(SAC_p0_dict)
u0 = extract_u0(SAC_u0_dict)
prob_func(du, u, p, t) = SAC_ODE_INH_EXC(du, u, p, t)
prob = SDEProblem(prob_func, noise1D, u0, tspan, p0)

n_traces = 20
initial_conditions = LinRange(-50, 40, n_traces)
function prob_func(prob, i, repeat; idx = 2)
     pI = prob.p
     pI[idx] = initial_conditions[i]
     println(pI)
     remake(prob, p = pI)
end

ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
@time sim = solve(ensemble_prob, SOSRI(), EnsembleDistributed(), trajectories = n_traces, 
     progress = true, progress_steps = 1, reltol = 0.01, abstol = 0.01, maxiters = 1e7)

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
ax11 = Axis(fSDE[2,4], title = "I_ext (pA)")

for sol in sim
     println(sol.t[end])
     Time = LinRange(sol.t[1], sol.t[end], 2000)
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
     lines!(ax11, Time, map(t -> sol(t)[1], Time))
end
display(fSDE)

#%% Running a Voltage clamp experiment
tspan = (0.0, 1e3)
SAC_p0_dict["VC"] = -60.0
SAC_p0_dict["g_GABA"] = 0.0
SAC_p0_dict["g_ACh"] = 0.0
p0 = extract_p0(SAC_p0_dict)
u0 = extract_u0(SAC_u0_dict)
prob_func(du, u, p, t) = SAC_ODE_VC(du, u, p, t; stim_start = 100.0, stim_stop = 500.0)
prob = SDEProblem(prob_func, noise1D, u0, tspan, p0)

n_traces = 20
initial_conditions = LinRange(-50, 40, n_traces)
function prob_func(prob, i, repeat; idx = 2)
     pI = prob.p
     pI[idx] = initial_conditions[i]
     println(pI)
     remake(prob, p = pI)
end

ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
@time sim = solve(ensemble_prob, SOSRI(), EnsembleDistributed(), trajectories = n_traces, 
     progress = true, progress_steps = 1, reltol = 0.01, abstol = 0.01, maxiters = 1e7)

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
ax11 = Axis(fSDE[2,4], title = "I_ext (pA)")

for sol in sim
     println(sol.t[end])
     Time = LinRange(sol.t[1], sol.t[end], 2000)
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
     lines!(ax11, Time, map(t -> sol(t)[1], Time))
end
display(fSDE)
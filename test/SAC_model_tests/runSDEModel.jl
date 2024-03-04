using PhysiologyModeling
using Pkg; Pkg.activate("test")
using PhysiologyPlotting
using GLMakie

#%% Run the example after here

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

tspan = (0.0, 100.0)
#Extract and modify parameters
p0_dict = SAC_p0_dict()
p0_dict["g_GABA"] = 0.0
p0_dict["g_ACh"] = 0.0
p0_dict["g_W"] = 0.075
p0 = extract_p0(p0_dict)
#Extract and modify initial conditions
u0_dict = SAC_u0_dict()
u0 = extract_u0(u0_dict)
#Set up the problem
prob_func(du, u, p, t) = SAC_ODE(du, u, p, t)
prob = SDEProblem(prob_func, noise1D, u0, tspan, p0)
@time sol = solve(prob, SOSRI(), reltol = 2e-2, abstol = 2e-2, progress = true, progress_steps = 1)

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
lines!(ax11, Time, map(t -> sol(t)[1], Time))


display(fSDE)
#%%
save("test/SAC_model_tests/data/SDESol.png", fSDE)


#%% Run an ensemble solution
tspan = (0.0, 60e3)
p0_dict = SAC_p0_dict()
#p0_dict["I_app"] = 10.0
p0_dict["g_K"] = 4.4
p0_dict["g_GABA"] = 0.0
p0_dict["g_ACh"] = 0.0

p0 = extract_p0(p0_dict)
u0_dict = SAC_u0_dict()
u0 = extract_u0(u0_dict)
prob_func(du, u, p, t) = SAC_ODE(du, u, p, t)
prob = SDEProblem(prob_func, noise1D, u0, tspan, p0)

n_traces = 30
initial_conditions = LinRange(4.0, 4.5, n_traces)
function prob_func(prob, i, repeat; idx = 4)
     pI = prob.p
     pI[idx] = initial_conditions[i]
     println(pI)
     remake(prob, p = pI) 
end

idx = findfirst(keys_p0 .== "g_Ca")
ensemble_prob = EnsembleProblem(prob, prob_func = (prob, i, repeat) -> prob_func(prob, i, repeat; idx = idx))
@time sim = solve(ensemble_prob, SOSRI(), EnsembleDistributed(), trajectories = n_traces, 
     progress = true, progress_steps = 1, reltol = 0.01, abstol = 0.01, maxiters = 1e7)

#
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
Colorbar(fSDE[5,1], limits = (initial_conditions[1], initial_conditions[end]), colormap = :viridis,
    vertical = false)
for (i, sol) in enumerate(sim)
     println(sol.t[end])
     Time = LinRange(sol.t[1], sol.t[end], 2000)
     x = initial_conditions[i]
     ln1 = lines!(ax1, Time, map(t -> sol(t)[2], Time), color = x, colormap = :viridis, colorrange = (initial_conditions[1], initial_conditions[end]))
     lines!(ax2, Time, map(t -> sol(t)[3], Time), color = x, colormap = :viridis, colorrange = (initial_conditions[1], initial_conditions[end]))
     lines!(ax3, Time, map(t -> sol(t)[4], Time), color = x, colormap = :viridis, colorrange = (initial_conditions[1], initial_conditions[end]))
     lines!(ax4, Time, map(t -> sol(t)[5], Time), color = x, colormap = :viridis, colorrange = (initial_conditions[1], initial_conditions[end]))
     lines!(ax5, Time, map(t -> sol(t)[6], Time), color = x, colormap = :viridis, colorrange = (initial_conditions[1], initial_conditions[end]))
     lines!(ax6, Time, map(t -> sol(t)[7], Time), color = x, colormap = :viridis, colorrange = (initial_conditions[1], initial_conditions[end]))
     lines!(ax7, Time, map(t -> sol(t)[8], Time), color = x, colormap = :viridis, colorrange = (initial_conditions[1], initial_conditions[end]))
     lines!(ax8, Time, map(t -> sol(t)[9], Time), color = x, colormap = :viridis, colorrange = (initial_conditions[1], initial_conditions[end]))
     lines!(ax9, Time, map(t -> sol(t)[10], Time), color = x, colormap = :viridis, colorrange = (initial_conditions[1], initial_conditions[end]))
     lines!(ax10, Time, map(t -> sol(t)[11], Time), color = x, colormap = :viridis, colorrange = (initial_conditions[1], initial_conditions[end]))
     lines!(ax11, Time, map(t -> sol(t)[1], Time), color = x, colormap = :viridis, colorrange = (initial_conditions[1], initial_conditions[end]))
end

display(fSDE)
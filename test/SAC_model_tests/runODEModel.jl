using Revise
using PhysiologyModeling
using PhysiologyPlotting
using GLMakie

#%% Step 1. Set up all parameters for the ODE
tspan = (0.0, 300e3)

SAC_p0_dict["I_app"] = 10.0
SAC_p0_dict["g_ACh"] = 0.0
SAC_p0_dict["g_GABA"] = 0.0
SAC_p0_dict["g_K"] = 10.0

prob = ODEProblem(SAC_ODE, vals_u0, tspan, extract_p0(SAC_p0_dict))
@time sol = solve(prob, progress = true, progress_steps = 1)

time = sol.t
fODE = Figure(resolution = (1800, 800))
ax1 = Axis(fODE[1,1], title = "Voltage (Vt)")
ax2 = Axis(fODE[2,1], title = "K Repol. (Nt)")
ax3 = Axis(fODE[3,1], title = "Na Gating (Mt)")
ax4 = Axis(fODE[4,1], title = "Na Close (Ht)")

ax5 = Axis(fODE[1,2], title = "Calcium (Ct)")
ax6 = Axis(fODE[2,2], title = "cAMP (At)")
ax7 = Axis(fODE[3,2], title = "TREK (Bt)")

ax8 = Axis(fODE[1,3], title = "ACh (Et)")
ax9 = Axis(fODE[2,3], title = "GABA (It)")

ax10 = Axis(fODE[1,4], title = "Noise (Wt)")

lines!(ax1, time, map(t -> sol(t)[1], time))
lines!(ax2, time, map(t -> sol(t)[2], time))
lines!(ax3, time, map(t -> sol(t)[3], time))
lines!(ax4, time, map(t -> sol(t)[4], time))
lines!(ax5, time, map(t -> sol(t)[5], time))
lines!(ax6, time, map(t -> sol(t)[6], time))
lines!(ax7, time, map(t -> sol(t)[7], time))
lines!(ax8, time, map(t -> sol(t)[8], time))
lines!(ax9, time, map(t -> sol(t)[9], time))
lines!(ax10, time, map(t -> sol(t)[10], time))
display(fODE)
save("test/SAC_model_tests/ODESol.png", fODE)
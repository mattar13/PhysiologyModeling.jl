using Revise, BenchmarkTools
using PhysiologyModeling
using PhysiologyPlotting
using GLMakie
using SparseArrays

#We want to make a map of SAC compartments
#%% Step 1. Set up all parameters for the ODE
tspan = (0.0, 10e3)

SAC_p0_dict["I_app"] = 0.0
SAC_p0_dict["g_ACh"] = 0.0
SAC_p0_dict["g_GABA"] = 0.0
SAC_p0_dict["g_K"] = 10.0

# Make the new model
n_compartments = 25
u0 = vcat(fill(vals_u0, n_compartments)'...)
SAC_p0_dict["C_m"] = 13.6#/n_compartments #Divide the cells in 5 parts

gK = LinRange(10.0, 1, n_compartments)
stim_starts = LinRange(0, 1000.0, n_compartments)
stim_dur = 500
fSAC(du, u, p, t) = SAC_ODE_Compartment(du, u, p, t; gK = gK, stim_starts = stim_starts, stim_dur = stim_dur)
prob = SDEProblem(fSAC, noise2D, u0, tspan, extract_p0(SAC_p0_dict))
@time sol = solve(prob, progress = true, progress_steps = 1)

time = sol.t[1]:10.0:sol.t[end]

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

for i in [1,n_compartments]
     println(i)
     lines!(ax1, time, map(t -> sol(t)[i, 1], time))
     lines!(ax2, time, map(t -> sol(t)[i, 2], time))
     lines!(ax3, time, map(t -> sol(t)[i, 3], time))
     lines!(ax4, time, map(t -> sol(t)[i, 4], time))
     lines!(ax5, time, map(t -> sol(t)[i, 5], time))
     lines!(ax6, time, map(t -> sol(t)[i, 6], time))
     lines!(ax7, time, map(t -> sol(t)[i, 7], time))
     lines!(ax8, time, map(t -> sol(t)[i, 8], time))
     lines!(ax9, time, map(t -> sol(t)[i, 9], time))
     lines!(ax10, time, map(t -> sol(t)[i, 10], time))
end
display(fODE)
save("test/SAC_model_tests/CompartmentSol_REV.png", fODE)
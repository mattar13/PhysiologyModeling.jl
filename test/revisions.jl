using Revise
using DifferentialEquations
using PhysiologyModeling
using PhysiologyPlotting
using GLMakie

#%%  Set up the cell map
domain_x = (xmin, xmax) = (0.0, 0.5)
domain_y = (ymin, ymax) = (0.0, 0.5)
dx = dy = 0.05 #Mean distribution is 40-50 micron (WR taylor et al)

cells = even_map(xmin = xmin, dx = dx, xmax = xmax, ymin = ymin, dy = dy, ymax = ymax)
#radii = rand(0.10:0.01:0.20, size(cells, 1))
radii = fill(0.200, size(cells, 1))
cell_map = CellMap(cells, radii, max_strength = 0.005);

# Step 1. Set up all parameters for the ODE
tspan = (0.0, 10e3)

SAC_p0_dict["I_app"] = 0.0
SAC_p0_dict["g_GABA"] = 0.0
SAC_p0_dict["g_ACh"] = 0.0
SAC_p0_dict["g_K"] = GK

p0 = extract_p0(SAC_p0_dict)
prob = SDEProblem(SAC_ODE_STIM, noise1D, vals_u0, tspan, p0)
@time sol = solve(prob, SOSRI(), progress = true, progress_steps = 1)

#%%
fSDE = Figure(resolution = (1800, 800))
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
display(fSDE)
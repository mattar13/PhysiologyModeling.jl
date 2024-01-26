using PhysiologyModeling

using Pkg; Pkg.activate("test")
using PhysiologyPlotting
using GLMakie

#%% Run the example after here
fSDE = Figure(size = (1000, 800))
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
ax11 = Axis(fSDE[2, 4])
# Step 1. Set up all parameters for the ODE
tspan = (0.0, 1000.0)

SAC_p0_dict["I_app"] = 10.0
SAC_p0_dict["g_GABA"] = 0.0
SAC_p0_dict["g_ACh"] = 0.0
SAC_p0_dict["g_TREK"] = 0.0
SAC_p0_dict["g_W"] = 0.0

param_modify = LinRange(2.0, 5.0, 20)
for (idx, x) in enumerate(param_modify)
     #SAC_p0_dict["I_app"] = x
     #SAC_p0_dict["g_W"] = x #may cause issues
     #SAC_p0_dict["g_TREK"] = x
     SAC_p0_dict["a_n"] = x
     p0 = extract_p0(SAC_p0_dict)
     prob = SDEProblem(SAC_ODE, noise1D, vals_u0, tspan, p0)
     #prob = ODEProblem(SAC_ODE, vals_u0, tspan, p0)
     @time sol = solve(prob, SOSRI(), reltol = 0.01, abstol = 0.01, progress = true, progress_steps = 1)
     #@time sol = solve(prob, reltol = 0.01, abstol = 0.01, progress = true, progress_steps = 1)
     # Plot the solution
     time = sol.t
     color_arr = [x]
     offset = 0.0
     println(color_arr)
     lines!(ax1, time, map(t -> sol(t)[1], time) .+ (idx * offset), colormap = :viridis, color = color_arr, colorrange = (param_modify[1], param_modify[end]))
     lines!(ax2, time, map(t -> sol(t)[2], time) .+ (idx * offset), colormap = :viridis, color = color_arr, colorrange = (param_modify[1], param_modify[end]))
     lines!(ax3, time, map(t -> sol(t)[3], time) .+ (idx * offset), colormap = :viridis, color = color_arr, colorrange = (param_modify[1], param_modify[end]))
     lines!(ax4, time, map(t -> sol(t)[4], time) .+ (idx * offset), colormap = :viridis, color = color_arr, colorrange = (param_modify[1], param_modify[end]))
     lines!(ax5, time, map(t -> sol(t)[5], time) .+ (idx * offset), colormap = :viridis, color = color_arr, colorrange = (param_modify[1], param_modify[end]))
     lines!(ax6, time, map(t -> sol(t)[6], time) .+ (idx * offset), colormap = :viridis, color = color_arr, colorrange = (param_modify[1], param_modify[end]))
     lines!(ax7, time, map(t -> sol(t)[7], time) .+ (idx * offset), colormap = :viridis, color = color_arr, colorrange = (param_modify[1], param_modify[end]))
     lines!(ax8, time, map(t -> sol(t)[8], time) .+ (idx * offset), colormap = :viridis, color = color_arr, colorrange = (param_modify[1], param_modify[end]))
     lines!(ax9, time, map(t -> sol(t)[9], time) .+ (idx * offset), colormap = :viridis, color = color_arr, colorrange = (param_modify[1], param_modify[end]))
     lines!(ax10, time, map(t -> sol(t)[10], time) .+ (idx * offset), colormap = :viridis, color = color_arr, colorrange = (param_modify[1], param_modify[end]))
     lines!(ax11, fill(x, length(time)), map(t -> sol(t)[1], time) .+ (idx * offset), colormap = :viridis, color = color_arr, colorrange = (param_modify[1], param_modify[end]))

     display(fSDE)
end
#save("test/SAC_model_tests/SDESol.png", fSDE)
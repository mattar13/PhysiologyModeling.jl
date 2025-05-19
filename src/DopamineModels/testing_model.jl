using Pkg; Pkg.activate(".")
using PhysiologyModeling
using Pkg; Pkg.activate("test")
using GLMakie

#%% ---------- load packages ----------
include("parameters.jl")
include("auxillary_functions.jl")
include("models.jl")

#%% ---------- simulation ----------
u0 = [-70.0, 0.05, 0.0, 0.0, 0.0, 2.0, 0.2]   # initial conditions
tspan = (0.0, 500.0)                          # ms
prob  = ODEProblem(dopaminergic_autoreceptor!, u0, tspan, p)
sol   = solve(prob, Tsit5(); saveat=0.5)      # 0.5 ms sampling

#%% ---------- plot ----------
fig = Figure()
ax1 = Axis(fig[1, 1], title="Dopaminergic Autoreceptor Model", xlabel="Time (ms)", ylabel="Membrane Potential (mV)")
ax2 = Axis(fig[1, 2], title="Calcium Concentration", xlabel="Time (ms)", ylabel="Calcium (mM)")
ax3 = Axis(fig[2, 1], title="Dopamine Concentration", xlabel="Time (ms)", ylabel="Dopamine (mM)")
ax4 = Axis(fig[2, 2], title="D2 Occupancy", xlabel="Time (ms)", ylabel="Occupancy")
ax5 = Axis(fig[3, 1], title="Gi/Go", xlabel="Time (ms)", ylabel="Gi/Go")
ax6 = Axis(fig[3, 2], title="cAMP", xlabel="Time (ms)", ylabel="cAMP")
ax7 = Axis(fig[4, 1], title="PKA", xlabel="Time (ms)", ylabel="PKA")

lines!(ax1, sol.t, sol[1, :], color=:blue, label="Membrane Potential (mV)")
lines!(ax2, sol.t, sol[2, :], color=:red, label="Calcium Concentration (mM)")
lines!(ax3, sol.t, sol[3, :], color=:green, label="Dopamine (mM)")
lines!(ax4, sol.t, sol[4, :], color=:orange, label="D2 Occupancy")
lines!(ax5, sol.t, sol[5, :], color=:purple, label="Gi/Go")
lines!(ax6, sol.t, sol[6, :], color=:cyan, label="cAMP")
lines!(ax7, sol.t, sol[7, :], color=:magenta, label="PKA")

display(fig)



#%% Test named tuple and storing parameters

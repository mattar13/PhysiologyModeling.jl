#=[Import packages]============================================================#
using ElectroPhysiology, PhysiologyModeling
using PhysiologyPlotting, GLMakie
using DiffEqCallbacks

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
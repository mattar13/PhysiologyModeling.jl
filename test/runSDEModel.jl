#=[Import packages]============================================================#
using ElectroPhysiology, PhysiologyModeling
using PhysiologyPlotting, GLMakie
using DiffEqCallbacks

#%%=[Solving a single SDE for tspan]==========================================#
#Extract and modify initial conditions
u0_dict = SAC_u0_dict()
u0_dict["a"] = 1.0
u0_dict["b"] = 0.0
u0 = extract_u0(u0_dict)

#Specify the timespan
tspan = (0.0, 300e3)
dt = 1.0
tseries = tspan[1]:dt:tspan[2]

#Set the stimulus parameters
save_fn = "D:/Data/Analysis/2024_02_12_WT_Cell1_IC.png"
p0_dict = SAC_p0_dict() #Extract parameters
p0_dict["g_GABA"] = 0.0 #Null GABA conductance
p0_dict["g_ACh"] = 0.0 #Null ACh conductance
#p0_dict["g_TREK"] = 0.0
p0 = extract_p0(p0_dict)

#create a callback for the glutamate pulse
function affect!(integrator; t1 = 0.0, dtpulse = 2500.0, tend = 10e3)
    dt = 0.0
    for t0 in t1:dtpulse:tend
        dt += gauss_pulse(integrator.t; t0 = t0, peak_amp = 0.05)
    end
    integrator.u[11] = dt
end
cb = PresetTimeCallback(tseries, affect!)

#Set up the problem
prob = SDEProblem(SAC_ODE, noise1D, u0, tspan, p0)
@time sol = solve(prob, SOSRI(), reltol = 2e-2, abstol = 2e-2, progress = true, progress_steps = 1)

# [Plot the solution]_________________________________________________________________________________________________________#
fig1 = Figure(size = (1800, 800))
ax1a = Axis(fig1[1,1], title = "I_ext (pA)")
ax2a = Axis(fig1[2,1], title = "Voltage (Vt)")
ax3a = Axis(fig1[3,1], title = "Noise (Wt)")

ax1b = Axis(fig1[1,2], title = "K Repol. (Nt)")
ax2b = Axis(fig1[2,2], title = "Na Gating (Mt)")
ax3b = Axis(fig1[3,2], title = "Na Close (Ht)")

ax1c = Axis(fig1[1,3], title = "Calcium (Ct)")
ax4d = Axis(fig1[2,3], title = "Dopamine (Dt)")
ax2c = Axis(fig1[3,3], title = "cAMP (At)")
ax3c = Axis(fig1[4,3], title = "TREK (Bt)")

ax1d = Axis(fig1[1,4], title = "ACh (Et)")
ax2d = Axis(fig1[2,4], title = "GABA (It)")
ax3d = Axis(fig1[3,4], title = "Glutamate (Gt)")

Time = sol.t
lines!(ax1a, Time, map(t -> sol(t)[1], Time))
lines!(ax2a, Time, map(t -> sol(t)[2], Time))
lines!(ax3a, Time, map(t -> sol(t)[13], Time))

lines!(ax1b, Time, map(t -> sol(t)[3], Time))
lines!(ax2b, Time, map(t -> sol(t)[4], Time))
lines!(ax3b, Time, map(t -> sol(t)[5], Time))

lines!(ax1c, Time, map(t -> sol(t)[6], Time))
lines!(ax2c, Time, map(t -> sol(t)[7], Time))
lines!(ax3c, Time, map(t -> sol(t)[8], Time))

lines!(ax1d, Time, map(t -> sol(t)[9], Time))
lines!(ax2d, Time, map(t -> sol(t)[10], Time))
lines!(ax3d, Time, map(t -> sol(t)[11], Time))
fig1
#%%
save("Modeling/Results/SDESol.png", fig1)
display(fig1)
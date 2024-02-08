using Revise
using ElectroPhysiology
using PhysiologyModeling
PhysiologyModeling.__init__()
using Pkg; Pkg.activate("test") #Activate the testing environment
using PhysiologyPlotting, GLMakie

#%%Run the ExcInc model according to this funciton
tspan = (0.0, 1000.0)
dt = 0.001
tstops = tspan[1]:dt:tspan[end]
input_wave = map(t -> gaussian(t), tstops)
InhExc_p0_dict["g_Exc"] = 1.0 #Set this experimentally
InhExc_p0_dict["g_Inh"] = 1.0 #Set this experimentally
p0 = extract_dict(InhExc_p0_dict, InhExc_p0_keys)
u0 = [0.0, 0.0]
prob = ODEProblem(InhExcModel, u0, tspan, p0)
@time sol = solve(prob, tstops = tstops, progress = true, progress_steps = 1)

fInhExc = Figure(size = (400, 400))
ax1 = Axis(fInhExc[1,1])
ax2 = Axis(fInhExc[2,1])
Time = sol.t
lines!(ax1, Time, map(t -> sol(t)[1], Time))
lines!(ax2, Time, map(t -> sol(t)[2], Time))
display(fInhExc)
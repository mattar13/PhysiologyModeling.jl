using Revise
using ElectroPhysiology
using PhysiologyModeling
PhysiologyModeling.__init__()
using Pkg; Pkg.activate("test") #Activate the testing environment
using PhysiologyPlotting, GLMakie

#%%Run the ExcInc model according to this funciton
tspan = (0.0, 1000.0)
dt = 0.1
tstops = tspan[1]:dt:tspan[end]
gE_input = rand(size(tstops))
gI_input = rand(size(tstops))

using MAT

file = matopen(raw"test\SAC_model_tests\eWT_PD.mat")
datafile = read(file)

p0 = extract_dict(InhExc_p0_dict, InhExc_p0_keys)
u0 = [0.0, 0.0, 0.0]
prob = ODEProblem(InhExcModel, u0, tspan, p0)
@time sol = solve(prob, tstops = tstops, progress = true, progress_steps = 1)

fInhExc = Figure(size = (400, 400))
ax_v = Axis(fInhExc[1,1])
ax_gE = Axis(fInhExc[2,1])
ax_gI = Axis(fInhExc[3,1])
Time = sol.t
lines!(ax_v, Time, map(t -> sol(t)[1], Time))
lines!(ax_gE, Time, map(t -> sol(t)[2], Time))
lines!(ax_gI, Time, map(t -> sol(t)[2], Time))

display(fInhExc)
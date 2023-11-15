using Revise
using BenchmarkTools, Profile, ProfileSVG
using PhysiologyModeling
using PhysiologyPlotting
using GLMakie

#%% First we want to benckmark the basic function 
du = similar(vals_u0)
u = vals_u0
p = extract_p0(SAC_p0_dict)
@benchmark SAC_ODE(du, u, p, 0.0)
@profile SAC_ODE(du, u, p, 0.0)
Profile.print()

#%% Profile the ODE model
tspan = (0.0, 300e3)

SAC_p0_dict["I_app"] = 10.0
SAC_p0_dict["g_ACh"] = 0.0
SAC_p0_dict["g_GABA"] = 0.0

prob = ODEProblem(SAC_ODE, vals_u0, tspan, extract_p0(SAC_p0_dict))
@time sol = solve(prob, progress = true, progress_steps = 1)
@benchmark sol = solve(prob, progress = true, progress_steps = 1)
@profile sol = solve(prob, progress = true, progress_steps = 1)
Profile.print()
ProfileSVG.save("ODE_Profile.svg")
Profile.clear()
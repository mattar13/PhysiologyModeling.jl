using Revise
using PhysiologyModeling
using PhysiologyPlotting
using GLMakie
using SparseArrays
import PhysiologyModeling.Φe

#%%
f_VS = Figure()
ax_vs = Axis(f_VS[1,1])

vs = -80:1:10
ϕe = map(v -> Φe(v, SAC_p0_dict["VSe"], SAC_p0_dict["V0e"])*SAC_p0_dict["ρe"], vs)
ϕi = map(v -> Φe(v, SAC_p0_dict["VSi"], SAC_p0_dict["V0i"])*SAC_p0_dict["ρi"], vs)

lines!(ax_vs, vs, ϕe)
lines!(ax_vs, vs, ϕi)

display(f_VS)
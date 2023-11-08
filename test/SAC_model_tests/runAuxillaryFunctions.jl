using Revise
using PhysiologyModeling
using PhysiologyPlotting
using GLMakie
using SparseArrays
import PhysiologyModeling: Φe, δX, IACh

#%% Plot out the acetylcholine release functions
f_REL = Figure()
ax_REL = Axis(f_REL[1,1], xlabel = "Voltage (mV)", ylabel = "NT Release (uM)")

VSe = SAC_p0_dict["VSe"]
VSi = SAC_p0_dict["VSi"]
V0e = SAC_p0_dict["V0e"]
V0i = SAC_p0_dict["V0i"]
ρe = SAC_p0_dict["ρe"]
ρi = SAC_p0_dict["ρi"]

vs = -80:1:10
ϕe = map(v -> Φe(v, VSe, V0e)*ρe, vs)
ϕi = map(v -> Φe(v, VSi, V0i)*ρi,  vs)

lines!(ax_REL, vs, ϕe)
lines!(ax_REL, vs, ϕi)

display(f_REL)

#%% Plot out the diffusion distance function
f_DIST = Figure()
ax_DIST = Axis(f_DIST[1,1], xlabel = "Distance from Soma (um)", ylabel = "NT Release (uM)")

max_strength = 0.05
max_dist = 0.15
slope_strength = 0.01

xs = 0.0:0.0005:0.2
δe = map(x -> δX(x, max_strength, max_dist, slope_strength), xs)

lines!(ax_DIST, xs, δe)
display(f_DIST)

#%% Plot out the current induced by neurotransmitters. The cell is clamped at
f_IVe = Figure()
ax_IVe = Axis(f_IVe[1,1], xlabel = "Voltage (mV)", ylabel = "Current (pA)")

vs = -80:1:10
es = 0.0:0.1:6.0
g_ACh = SAC_p0_dict["g_ACh"]
k_ACh = SAC_p0_dict["k_ACh"]
E_ACh = SAC_p0_dict["E_ACh"]

for e in es
     Ie = map(v -> IACh(v, e, g_ACh, k_ACh, E_ACh), vs)
     lines!(ax_IVe, vs, Ie)
end
display(f_IVe)
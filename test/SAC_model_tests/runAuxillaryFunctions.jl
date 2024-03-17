using Revise
using PhysiologyModeling
import PhysiologyModeling: ring, Φe, IACh, IGABA, ħe, ħi, ring_circle_overlap_area
using Pkg; Pkg.activate("test")

using PhysiologyPlotting
using GLMakie
using SparseArrays


#%% 
f = Figure()
ax1 = Axis(f[1,1])
tseries = 0.0:0.1:100.0
lines!(ax1, tseries, gauss_pulse.(tseries, t0 = 10.0))
display(f)


#%%
f_DIST = Figure()
ax_DIST = Axis(f_DIST[1,1], xlabel = "Distance from Soma (um)", ylabel = "NT Release (uM)")

xs = 0.0:0.001:0.30
for rIN in LinRange(0.01, 0.1, 30)
     dist_func1(d) = ring_circle_overlap_area(d; density = 0.01, r_inner = rIN, r_outer = 0.18, r_circle = 0.18)
     dist_func2(d) = ring(d; max_strength = 0.001, max_dist = 0.18, slope = 0.025)
     δe1 = map(d -> dist_func1(d), xs)
     δe2 = map(d -> dist_func2(d), xs)

     lines!(ax_DIST, xs, δe1)
     lines!(ax_DIST, xs, δe2)
end

display(f_DIST)
save("test/SAC_model_tests/data/DistanceFunc.png", f_DIST)

#%% Plot out the current induced by neurotransmitters. The cell is clamped at
f_REL = Figure(resolution = (1200, 800))

ax_Ne = Axis(f_REL[1,1], title = "Voltage gated NT Rel", xlabel = "Voltage (mV)", ylabel = "NT Release (uM)")
ax_Ni = Axis(f_REL[2,1], xlabel = "Voltage (mV)", ylabel = "NT Release (uM)")

ax_NHe = Axis(f_REL[1,2], title = "NT channel activation", xlabel = "neurotransmitters (μM)", ylabel = "Norm. Act.")
ax_NHi = Axis(f_REL[2,2], xlabel = "neurotransmitters (μM)", ylabel = "Norm. Act.")

ax_IVe = Axis(f_REL[1,3], title = "Current induced by NT", xlabel = "Voltage (mV)", ylabel = "Current (pA)")
ax_IVi = Axis(f_REL[2,3], xlabel = "Voltage (mV)", ylabel = "Current (pA)")

vs = -80:1:10
es = 0.0:0.05:5.0
is = 0.0:0.05:5.0

VSe = SAC_p0_dict["VSe"]
VSi = SAC_p0_dict["VSi"]
V0e = SAC_p0_dict["V0e"]
V0i = SAC_p0_dict["V0i"]
ρe = SAC_p0_dict["ρe"]
ρi = SAC_p0_dict["ρi"]

g_ACh = SAC_p0_dict["g_ACh"]
k_ACh = SAC_p0_dict["k_ACh"]
E_ACh = SAC_p0_dict["E_ACh"]

g_GABA = SAC_p0_dict["g_GABA"]
k_GABA = SAC_p0_dict["k_GABA"]
E_Cl = SAC_p0_dict["E_Cl"]

ϕe = map(v -> Φe(v, VSe, V0e)*ρe, vs)
ϕi = map(v -> Φe(v, VSi, V0i)*ρi,  vs)

lines!(ax_Ne, vs, ϕe, colormap = :thermal, color = ϕe, colorrange = (0.0, maximum(es)))
lines!(ax_Ni, vs, ϕi, colormap = :thermal, color = ϕi, colorrange = (0.0, maximum(es)))

He = map(e-> ħe(e, k_ACh), es)
Hi = map(i-> ħe(i, k_GABA), is)

lines!(ax_NHe, es, He, colormap = :thermal, color = es, colorrange = (0.0, maximum(es)))
lines!(ax_NHi, is, Hi, colormap = :thermal, color = is, colorrange = (0.0, maximum(is)))

for e in es
     Ie = map(v -> IACh(v, e, g_ACh, k_ACh, E_ACh), vs)
     lines!(ax_IVe, vs, Ie, colormap = :thermal, color = [e], colorrange = (0.0, maximum(es)))
end

for i in is
     Ii = map(v -> IGABA(v, i, g_GABA, k_GABA, E_Cl), vs)
     lines!(ax_IVi, vs, Ii, colormap = :thermal, color = [i], colorrange = (0.0, maximum(is)))
end
save("test/SAC_model_tests/data/Diffusion.png", f_REL)
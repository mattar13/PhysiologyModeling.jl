using Pkg; Pkg.activate("PhysiologyModeling.jl")
using Catalyst
using OrdinaryDiffEq, SteadyStateDiffEq
using GLMakie
using ModelingToolkit

include("parameters.jl")
include("auxillary_functions.jl")
include("models.jl")

# Time span
tspan = (0.0, 5000.0)
# Create and solve the problem

ss_prob = SteadyStateProblem(metabolic_network, u0, ss_params)
ss_sol = solve(ss_prob, DynamicSS(Rodas5()))
u0 = ss_sol.u

prob = ODEProblem(metabolic_network, u0, tspan, params)
# Use a more robust solver for stiff systems
sol = solve(prob, Rodas5(), abstol=1e-8, reltol=1e-8, maxiters=100000)

for i in 1:length(metabolic_network.unknowns)
    println("Index $i: ", metabolic_network.unknowns[i])
end
# Extract time series for plotting
tseries = LinRange(tspan[1], tspan[2], 10000)

# Extract time series using correct variable order
# Index mapping based on actual solver output:
# 1: ATP, 2: ADP, 3: Pi, 4: AMP, 5: Ado, 6: GLU, 7: NAD, 8: Pyruvate, 9: NADH, 10: CoA, 11: Acetyl_CoA, 12: CO2, 13: FAD, 14: FADH2, 15: O2, 16: H2O, 17: Lactate, 18: Alanine, 19: Cr, 20: PCr, 21: H2CO3, 22: HCO3, 23: H, 24: V, 25: h, 26: m, 27: n, 28: a, 29: d
atp_series        = map(t -> sol(t)[1], tseries)
adp_series        = map(t -> sol(t)[2], tseries)
pi_series         = map(t -> sol(t)[3], tseries)
amp_series        = map(t -> sol(t)[4], tseries)
ado_series        = map(t -> sol(t)[5], tseries)
glucose_series    = map(t -> sol(t)[6], tseries)
nad_series        = map(t -> sol(t)[7], tseries)
pyruvate_series   = map(t -> sol(t)[8], tseries)
nadh_series       = map(t -> sol(t)[9], tseries)
coa_series        = map(t -> sol(t)[10], tseries)
acetyl_coa_series = map(t -> sol(t)[11], tseries)
co2_series        = map(t -> sol(t)[12], tseries)
fad_series        = map(t -> sol(t)[13], tseries)
fadh2_series      = map(t -> sol(t)[14], tseries)
oxygen_series     = map(t -> sol(t)[15], tseries)
h2o_series        = map(t -> sol(t)[16], tseries)
lactate_series    = map(t -> sol(t)[17], tseries)
alanine_series    = map(t -> sol(t)[18], tseries)
cr_series         = map(t -> sol(t)[19], tseries)
pcr_series        = map(t -> sol(t)[20], tseries)
h2co3_series      = map(t -> sol(t)[21], tseries)
hco3_series       = map(t -> sol(t)[22], tseries)
h_ion_series      = map(t -> sol(t)[23], tseries)
v_series          = map(t -> sol(t)[24], tseries)
h_series          = map(t -> sol(t)[25], tseries)
m_series          = map(t -> sol(t)[26], tseries)
n_series          = map(t -> sol(t)[27], tseries)
a_series          = map(t -> sol(t)[28], tseries)
d_series          = map(t -> sol(t)[29], tseries)

# Plotting
fig = Figure(size = (1200, 1000))

# Row 1: Nucleotides and energy
ax1a = Axis(fig[1, 1], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "ATP", ytickformat = "{:.2f}")
ax1b = Axis(fig[1, 2], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "ADP", ytickformat = "{:.2f}")
ax1c = Axis(fig[1, 3], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "AMP", ytickformat = "{:.2f}")
ax1d = Axis(fig[1, 4], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "Ado", ytickformat = "{:.2f}")
ax1e = Axis(fig[1, 5], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "Pi", ytickformat = "{:.2f}")

# Row 2: R1P and PCr/Cr
ax2a = Axis(fig[2, 1], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "PCr & Cr", ytickformat = "{:.2f}")
ax2b = Axis(fig[2, 2], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "PCr & Cr", ytickformat = "{:.2f}")
# Row 3: GLU and metabolic intermediates
ax3a = Axis(fig[3, 1], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "GLU", ytickformat = "{:.2f}")
ax3b = Axis(fig[3, 2], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "Pyruvate", ytickformat = "{:.2f}")
ax3c = Axis(fig[3, 3], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "Lactate", ytickformat = "{:.2f}")
ax3d = Axis(fig[3, 4], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "Alanine", ytickformat = "{:.2f}")

# Row 4: TCA cycle intermediates
ax4a = Axis(fig[4, 1], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "Acetyl-CoA", ytickformat = "{:.2f}")
ax4b = Axis(fig[4, 2], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "CoA", ytickformat = "{:.2f}")

# Row 5: Metabolic and redox compounds
ax5a = Axis(fig[5, 1], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "O2, CO2 & H2O", ytickformat = "{:.2f}")
ax5b = Axis(fig[5, 2], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "FAD & FADH2", ytickformat = "{:.2f}")
ax5c = Axis(fig[5, 3], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "NAD & NADH", ytickformat = "{:.2f}")

# Row 6: Carbonic acid system
ax6a = Axis(fig[6, 1], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "H2CO3", ytickformat = "{:.2f}")
ax6b = Axis(fig[6, 2], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "HCO3-", ytickformat = "{:.2f}")
ax6c = Axis(fig[6, 3], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "H+ (×1000)", ytickformat = "{:.2f}")

# Row 7: Voltage and gating variables
ax7a = Axis(fig[7, 1], xlabel = "Time (ms)", ylabel = "Voltage (mV)", title = "Membrane Potential", ytickformat = "{:.2f}")
ax7b = Axis(fig[7, 2], xlabel = "Time (ms)", ylabel = "Gating Variable", title = "Na Activation (m)", ytickformat = "{:.2f}")
ax7c = Axis(fig[7, 3], xlabel = "Time (ms)", ylabel = "Gating Variable", title = "Na Inactivation (h)", ytickformat = "{:.2f}")
ax7d = Axis(fig[7, 4], xlabel = "Time (ms)", ylabel = "Gating Variable", title = "K Activation (n)", ytickformat = "{:.2f}")
ax7e = Axis(fig[7, 5], xlabel = "Time (ms)", ylabel = "Pump Activation (a)", title = "NaK Pump Activation", ytickformat = "{:.2f}")
ax7f = Axis(fig[7, 6], xlabel = "Time (ms)", ylabel = "Pump Activation (d)", title = "NaK Pump Drive", ytickformat = "{:.2f}")

# Plot all series
lines!(ax1a, tseries, atp_series, color = :blue, linewidth = 2)
lines!(ax1b, tseries, adp_series, color = :orange, linewidth = 2)
lines!(ax1c, tseries, amp_series, color = :green, linewidth = 2)
lines!(ax1d, tseries, ado_series, color = :red, linewidth = 2)
lines!(ax1e, tseries, pi_series, color = :purple, linewidth = 2)

lines_PCR = lines!(ax2a, tseries, pcr_series, color = :darkblue, linewidth = 2, label = "PCr")
lines_CR = lines!(ax2b, tseries, cr_series, color = :darkorange, linewidth = 2, label = "Cr")
legend_pcr_cr = Legend(fig[2, 3], [lines_PCR, lines_CR], ["PCr", "Cr"], tellwidth = false, tellheight = false)

lines!(ax3a, tseries, glucose_series, color = :brown, linewidth = 2)
lines!(ax3b, tseries, pyruvate_series, color = :orange, linewidth = 2)
lines!(ax3c, tseries, lactate_series, color = :pink, linewidth = 2)
lines!(ax3d, tseries, alanine_series, color = :lightblue, linewidth = 2)

lines!(ax4a, tseries, acetyl_coa_series, color = :olive, linewidth = 2)
lines!(ax4b, tseries, coa_series, color = :teal, linewidth = 2)

line_o2 = lines!(ax5a, tseries, oxygen_series, color = :cyan, linewidth = 2, label = "O2")
line_co2 = lines!(ax5a, tseries, co2_series, color = :gray, linewidth = 2, label = "CO2")
line_h2o = lines!(ax5a, tseries, h2o_series, color = :lightblue, linewidth = 2, label = "H2O")

line_fad = lines!(ax5b, tseries, fad_series, color = :khaki, linewidth = 2, label = "FAD")
line_fadh2 = lines!(ax5b, tseries, fadh2_series, color = :lavender, linewidth = 2, label = "FADH2")

line_nad = lines!(ax5c, tseries, nad_series, color = :forestgreen, linewidth = 2, label = "NAD")
line_nadh = lines!(ax5c, tseries, nadh_series, color = :hotpink, linewidth = 2, label = "NADH")

# Remove grid-based legends and use floating legends anchored to the right of each axis
legend_o2 = Legend(fig[5, 4], [line_o2, line_co2, line_h2o], ["O2", "CO2", "H2O"], tellwidth = false, tellheight = false)
legend_fad = Legend(fig[5, 5], [line_fad, line_fadh2, line_nad, line_nadh], ["FAD", "FADH2", "NAD", "NADH"], tellwidth = false, tellheight = false)

# Plot carbonic acid system
lines!(ax6a, tseries, h2co3_series, color = :darkblue, linewidth = 2)
lines!(ax6b, tseries, hco3_series, color = :blue, linewidth = 2)
lines!(ax6c, tseries, h_ion_series .* 1000, color = :red, linewidth = 2)  # Scale H+ by 1000 for visibility

# Plot voltage and gating variables
lines!(ax7a, tseries, v_series, color = :black, linewidth = 3)
lines!(ax7b, tseries, m_series, color = :red, linewidth = 2)
lines!(ax7c, tseries, h_series, color = :blue, linewidth = 2)
lines!(ax7d, tseries, n_series, color = :green, linewidth = 2)
lines!(ax7e, tseries, a_series, color = :purple, linewidth = 2)
lines!(ax7f, tseries, d_series, color = :orange, linewidth = 2)
fig
#%% Figure 2 Just the voltage 
fig2 = Figure(size = (800, 500))

ax1 = Axis(fig2[1, 1:4], xlabel = "Time (ms)", ylabel = "Voltage (mV)", title = "Membrane Potential", ytickformat = "{:.2f}")
ax2 = Axis(fig2[2, 1], xlabel = "Time (ms)", ylabel = "Gating Variable", title = "Na Activation (m)", ytickformat = "{:.2f}")
ax3 = Axis(fig2[2, 2], xlabel = "Time (ms)", ylabel = "Gating Variable", title = "Na Inactivation (h)", ytickformat = "{:.2f}")
ax4 = Axis(fig2[2, 3], xlabel = "Time (ms)", ylabel = "Gating Variable", title = "K Activation (n)", ytickformat = "{:.2f}")
ax5 = Axis(fig2[2, 4], xlabel = "Time (ms)", ylabel = "Pump Activation", title = "Pump Activation (a)", ytickformat = "{:.2f}")
lines!(ax1, tseries, v_series, color = :black, linewidth = 3)
lines!(ax2, tseries, m_series, color = :red, linewidth = 2)
lines!(ax3, tseries, h_series, color = :blue, linewidth = 2)
lines!(ax4, tseries, n_series, color = :green, linewidth = 2)
lines!(ax5, tseries, a_series, color = :purple, linewidth = 2)

#%% Plot the voltage vs the ATP
fig3 = Figure(size = (800, 500))

ax1 = Axis(fig3[1, 1], xlabel = "Time (ms)", ylabel = "Voltage (mV)", title = "Membrane Potential", ytickformat = "{:.2f}")
ax2 = Axis(fig3[2, 1], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "ATP", ytickformat = "{:.2f}")
ax3 = Axis(fig3[1:2, 2], xlabel = "Voltage (mV)", ylabel = "ATP (mM)", title = "ATP vs Voltage", xtickformat = "{:.2f}", ytickformat = "{:.2f}")

lines!(ax1, tseries, v_series, color = :black, linewidth = 3)
lines!(ax2, tseries, atp_series, color = :blue, linewidth = 2)
lines!(ax3, v_series, atp_series, color = :red, linewidth = 2)

fig3
#%% Plot the Pump activation
a_ATP_half_val = 0.005
a_V_half_val = -30.0
a_V_slope_val = 10.0

atp_rng = LinRange(0.0, 10.0, 1000)
v_rng = LinRange(-100.0, 0.0, 1000)
fATP = (atp_rng./(atp_rng .+ a_ATP_half_val))
fV = 1 ./ (1 .+ exp.(-(v_rng .- -30.0)/10.0))

fig4 = Figure(size = (800, 500))

ax1 = Axis(fig4[1, 1], xlabel = "ATP (mM)", ylabel = "Pump Activation (a)", title = "Pump Activation", xtickformat = "{:.2f}", ytickformat = "{:.2f}")
ax2 = Axis(fig4[1, 2], xlabel = "Voltage (mV)", ylabel = "Pump Activation (a)", title = "Pump Activation", xtickformat = "{:.2f}", ytickformat = "{:.2f}")

lines!(ax1, atp_rng, fATP, color = :black, linewidth = 3)
lines!(ax2, v_rng, fV, color = :blue, linewidth = 2)
fig4
#%% Which plot do you want to see?
fig
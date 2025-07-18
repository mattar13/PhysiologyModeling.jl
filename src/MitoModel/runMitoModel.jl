using Pkg; Pkg.activate("PhysiologyModeling.jl")
using Catalyst
using OrdinaryDiffEq
using GLMakie
using ModelingToolkit

include("parameters.jl")
include("auxillary_functions.jl")
include("models.jl")

# Time span
tspan = (0.0, 5000.0)
# Create and solve the problem
prob = ODEProblem(metabolic_network, u0, tspan, params)
# Use a more robust solver for stiff systems
sol = solve(prob, Rodas5(), abstol=1e-8, reltol=1e-8, maxiters=100000)

for i in 1:length(metabolic_network.unknowns)
    println("Index $i: ", metabolic_network.unknowns[i])
end
# Extract time series using indices (we'll need to determine the correct order)
# For now, let's use the first few variables and see what we get
tseries = LinRange(tspan[1], tspan[2], 1000)

# Extract time series using correct variable order
# Index mapping based on actual solver output (PPP removed):
# 1: ATP, 2: ADP, 3: Pi, 4: AMP, 5: Ado, 6: GLU, 7: NAD, 8: Pyruvate, 9: NADH, 10: CoA, 11: Acetyl_CoA, 12: CO2, 13: FAD, 14: FADH2, 15: O2, 16: H2O, 17: PEP, 18: BPG, 19: P3G, 20: Lactate, 21: Alanine, 22: Cr, 23: PCr, 24: V, 25: h, 26: m, 27: n, 28: a, 29: d

atp_series        = round.(map(t -> sol(t)[1], tseries), digits=place_digits)
adp_series        = round.(map(t -> sol(t)[2], tseries), digits=place_digits)
pi_series         = round.(map(t -> sol(t)[3], tseries), digits=place_digits)
amp_series        = round.(map(t -> sol(t)[4], tseries), digits=place_digits)
ado_series        = round.(map(t -> sol(t)[5], tseries), digits=place_digits)
glucose_series    = round.(map(t -> sol(t)[6], tseries), digits=place_digits)
nad_series        = round.(map(t -> sol(t)[7], tseries), digits=place_digits)
pyruvate_series   = round.(map(t -> sol(t)[8], tseries), digits=place_digits)
nadh_series       = round.(map(t -> sol(t)[9], tseries), digits=place_digits)
coa_series        = round.(map(t -> sol(t)[10], tseries), digits=place_digits)
acetyl_coa_series = round.(map(t -> sol(t)[11], tseries), digits=place_digits)
co2_series        = round.(map(t -> sol(t)[12], tseries), digits=place_digits)
fad_series        = round.(map(t -> sol(t)[13], tseries), digits=place_digits)
fadh2_series      = round.(map(t -> sol(t)[14], tseries), digits=place_digits)
oxygen_series     = round.(map(t -> sol(t)[15], tseries), digits=place_digits)
h2o_series        = round.(map(t -> sol(t)[16], tseries), digits=place_digits)
pep_series        = round.(map(t -> sol(t)[17], tseries), digits=place_digits)
bpg_series        = round.(map(t -> sol(t)[18], tseries), digits=place_digits)
p3g_series        = round.(map(t -> sol(t)[19], tseries), digits=place_digits)
lactate_series    = round.(map(t -> sol(t)[20], tseries), digits=place_digits)
alanine_series    = round.(map(t -> sol(t)[21], tseries), digits=place_digits)
cr_series         = round.(map(t -> sol(t)[22], tseries), digits=place_digits)
pcr_series        = round.(map(t -> sol(t)[23], tseries), digits=place_digits)
v_series          = round.(map(t -> sol(t)[24], tseries), digits=place_digits)
h_series          = round.(map(t -> sol(t)[25], tseries), digits=place_digits)
m_series          = round.(map(t -> sol(t)[26], tseries), digits=place_digits)
n_series          = round.(map(t -> sol(t)[27], tseries), digits=place_digits)
a_series          = round.(map(t -> sol(t)[28], tseries), digits=place_digits)
d_series          = round.(map(t -> sol(t)[29], tseries), digits=place_digits)

# Plotting
fig = Figure(size = (1200, 1000))

# Row 1: Nucleotides and energy
ax1a = Axis(fig[1, 1], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "ATP")
ax1b = Axis(fig[1, 2], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "ADP")
ax1c = Axis(fig[1, 3], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "AMP")
ax1d = Axis(fig[1, 4], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "Ado")
ax1e = Axis(fig[1, 5], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "Pi")

# Row 2: R1P and PCr/Cr
ax2a = Axis(fig[2, 1], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "PCr & Cr")
ax2b = Axis(fig[2, 2], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "PCr & Cr")
# Row 3: GLU and glycolysis elements
ax3a = Axis(fig[3, 1], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "GLU")
ax3b = Axis(fig[3, 2], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "PEP")
ax3c = Axis(fig[3, 3], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "BPG")
ax3d = Axis(fig[3, 4], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "P3G")
ax3e = Axis(fig[3, 5], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "Pyruvate")
ax3f = Axis(fig[3, 6], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "Lactate")
ax3g = Axis(fig[3, 7], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "Alanine")

# Row 4: TCA cycle intermediates
ax4a = Axis(fig[4, 1], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "3PG")
ax4b = Axis(fig[4, 2], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "Acetyl-CoA")
ax4c = Axis(fig[4, 3], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "CoA")

# Row 5: Metabolic and redox compounds
ax5a = Axis(fig[5, 1], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "O2, CO2 & H2O")
ax5b = Axis(fig[5, 2], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "FAD & FADH2")
ax5c = Axis(fig[5, 3], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "NAD & NADH")

# Row 6: Voltage and gating variables
ax6a = Axis(fig[6, 1], xlabel = "Time (ms)", ylabel = "Voltage (mV)", title = "Membrane Potential")
ax6b = Axis(fig[6, 2], xlabel = "Time (ms)", ylabel = "Gating Variable", title = "Na Activation (m)")
ax6c = Axis(fig[6, 3], xlabel = "Time (ms)", ylabel = "Gating Variable", title = "Na Inactivation (h)")
ax6d = Axis(fig[6, 4], xlabel = "Time (ms)", ylabel = "Gating Variable", title = "K Activation (n)")
ax6e = Axis(fig[6, 5], xlabel = "Time (ms)", ylabel = "Pump Activation (a)", title = "NaK Pump Activation")
ax6f = Axis(fig[6, 6], xlabel = "Time (ms)", ylabel = "Pump Activation (d)", title = "NaK Pump Drive")

# Plot all series
lines!(ax1a, tseries, atp_series, color = :blue, linewidth = 2)
lines!(ax1b, tseries, adp_series, color = :orange, linewidth = 2)
lines!(ax1c, tseries, amp_series, color = :green, linewidth = 2)
lines!(ax1d, tseries, ado_series, color = :red, linewidth = 2)
lines!(ax1e, tseries, pi_series, color = :purple, linewidth = 2)

lines!(ax2a, tseries, pcr_series, color = :darkblue, linewidth = 2, label = "PCr")
lines_CR = lines!(ax2b, tseries, cr_series, color = :darkorange, linewidth = 2, label = "Cr")
legend_pcr_cr = Legend(fig[2, 3], [lines_PCR, lines_CR], ["PCr", "Cr"], tellwidth = false, tellheight = false)

lines!(ax3a, tseries, glucose_series, color = :brown, linewidth = 2)
lines!(ax3b, tseries, pep_series, color = :darkgreen, linewidth = 2)
lines!(ax3c, tseries, bpg_series, color = :darkred, linewidth = 2)
lines!(ax3d, tseries, p3g_series, color = :brown, linewidth = 2)
lines!(ax3e, tseries, pyruvate_series, color = :orange, linewidth = 2)
lines!(ax3f, tseries, lactate_series, color = :pink, linewidth = 2)
lines!(ax3g, tseries, alanine_series, color = :lightblue, linewidth = 2)

lines!(ax4a, tseries, p3g_series, color = :brown, linewidth = 2)
lines!(ax4b, tseries, acetyl_coa_series, color = :olive, linewidth = 2)
lines!(ax4c, tseries, coa_series, color = :teal, linewidth = 2)

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


# Plot voltage and gating variables
lines!(ax6a, tseries, v_series, color = :black, linewidth = 3)
lines!(ax6b, tseries, m_series, color = :red, linewidth = 2)
lines!(ax6c, tseries, h_series, color = :blue, linewidth = 2)
lines!(ax6d, tseries, n_series, color = :green, linewidth = 2)
lines!(ax6e, tseries, a_series, color = :purple, linewidth = 2)
lines!(ax6f, tseries, d_series, color = :orange, linewidth = 2)
fig
#%% Figure 2 Just the voltage 
fig2 = Figure(size = (800, 500))

ax1 = Axis(fig2[1, 1:4], xlabel = "Time (ms)", ylabel = "Voltage (mV)", title = "Membrane Potential")
ax2 = Axis(fig2[2, 1], xlabel = "Time (ms)", ylabel = "Voltage (mV)", title = "Na Activation (m)")
ax3 = Axis(fig2[2, 2], xlabel = "Time (ms)", ylabel = "Voltage (mV)", title = "Na Inactivation (h)")
ax4 = Axis(fig2[2, 3], xlabel = "Time (ms)", ylabel = "Voltage (mV)", title = "K Activation (n)")
ax5 = Axis(fig2[2, 4], xlabel = "Time (ms)", ylabel = "Voltage (mV)", title = "Pump Activation (a)")
lines!(ax1, tseries, v_series, color = :black, linewidth = 3)
lines!(ax2, tseries, m_series, color = :red, linewidth = 2)
lines!(ax3, tseries, h_series, color = :blue, linewidth = 2)
lines!(ax4, tseries, n_series, color = :green, linewidth = 2)
lines!(ax5, tseries, a_series, color = :purple, linewidth = 2)

#%% Plot the voltage vs the ATP
fig3 = Figure(size = (800, 500))

ax1 = Axis(fig3[1, 1], xlabel = "Time (ms)", ylabel = "Voltage (mV)", title = "Membrane Potential")
ax2 = Axis(fig3[2, 1], xlabel = "Time (ms)", ylabel = "Voltage (mV)", title = "ATP")
ax3 = Axis(fig3[1:2, 2], xlabel = "Time (ms)", ylabel = "Voltage (mV)", title = "Cycle")

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

ax1 = Axis(fig4[1, 1], xlabel = "ATP (mM)", ylabel = "Pump Activation (a)", title = "Pump Activation")
ax2 = Axis(fig4[1, 2], xlabel = "Voltage (mV)", ylabel = "Pump Activation (a)", title = "Pump Activation")

lines!(ax1, atp_rng, fATP, color = :black, linewidth = 3)
lines!(ax2, v_rng, fV, color = :blue, linewidth = 2)
fig4
#%% Which plot do you want to see?
fig
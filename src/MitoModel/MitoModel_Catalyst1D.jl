using Pkg; Pkg.activate("PhysiologyModeling.jl")
using Catalyst
using OrdinaryDiffEq
using GLMakie
using ModelingToolkit

include("auxillary_functions.jl")
include("parameters.jl")
include("models.jl")

complete_system |> typeof |> fieldnames
complete_system.systems

# Time span
tspan = (0.0, 1000.0)
# Create and solve the problem
prob = ODEProblem(complete_system, u0, tspan, params)
# Use a more robust solver for stiff systems
sol = solve(prob, Rodas5(), abstol=1e-8, reltol=1e-8, maxiters=10000)

#Extract the variables
tseries = LinRange(tspan[1], tspan[2], 1000)

# Print solution structure to understand variable ordering
println("Solution structure:")
println("Number of variables: ", length(sol.u[1]))
println("First solution: ", sol.u[1])

# Extract time series using indices (we'll need to determine the correct order)
# For now, let's use the first few variables and see what we get
atp_series = map(t -> sol(t)[1], tseries)
adp_series = map(t -> sol(t)[2], tseries)
pi_series = map(t -> sol(t)[3], tseries)
amp_series = map(t -> sol(t)[4], tseries)
ado_series = map(t -> sol(t)[5], tseries)
ribose1p_series = map(t -> sol(t)[6], tseries)
glucose_series = map(t -> sol(t)[7], tseries)
oxygen_series = map(t -> sol(t)[8], tseries)
pcr_series = map(t -> sol(t)[9], tseries)
cr_series = map(t -> sol(t)[10], tseries)
co2_series = map(t -> sol(t)[11], tseries)
h2o_series = map(t -> sol(t)[12], tseries)
pep_series = map(t -> sol(t)[13], tseries)
bpg_series = map(t -> sol(t)[14], tseries)
pyruvate_series = map(t -> sol(t)[15], tseries)
p3g_series = map(t -> sol(t)[16], tseries)
acetyl_coa_series = map(t -> sol(t)[17], tseries)
coa_series = map(t -> sol(t)[18], tseries)
nad_series = map(t -> sol(t)[19], tseries)
nadh_series = map(t -> sol(t)[20], tseries)
fad_series = map(t -> sol(t)[21], tseries)
fadh2_series = map(t -> sol(t)[22], tseries)
v_series = map(t -> sol(t)[26], tseries)
m_series = map(t -> sol(t)[23], tseries)
h_series = map(t -> sol(t)[24], tseries)
n_series = map(t -> sol(t)[25], tseries)

# Plotting
fig = Figure(size = (1200, 1000))

# Row 1: Nucleotides and energy
ax1a = Axis(fig[1, 1], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "ATP")
ax1b = Axis(fig[1, 2], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "ADP")
ax1c = Axis(fig[1, 3], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "AMP")
ax1d = Axis(fig[1, 4], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "Ado")
ax1e = Axis(fig[1, 5], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "Pi")
ax1f = Axis(fig[1, 6], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "R1P")

# Row 2: Creatine system and substrates
ax2a = Axis(fig[2, 1], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "PCr & Cr")
ax2b = Axis(fig[2, 2], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "GLU")
ax2c = Axis(fig[2, 3], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "PEP")
ax2d = Axis(fig[2, 4], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "BPG")
ax2e = Axis(fig[2, 5], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "Pyruvate")

# Row 3: TCA cycle intermediates
ax3a = Axis(fig[3, 1], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "3PG")
ax3b = Axis(fig[3, 2], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "Acetyl-CoA")
ax3c = Axis(fig[3, 3], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "CoA")

# Row 4: Electron carriers
ax4a = Axis(fig[4, 1], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "NAD & NADH")
ax4b = Axis(fig[4, 2], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "FAD & FADH2")

# Row 5: Molecular compounds
ax5a = Axis(fig[5, 1], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "O2, CO2 & H2O")

# Row 6: Voltage and gating variables
ax6a = Axis(fig[6, 1], xlabel = "Time (ms)", ylabel = "Voltage (mV)", title = "Membrane Potential")
ax6b = Axis(fig[6, 2], xlabel = "Time (ms)", ylabel = "Gating Variable", title = "Na Activation (m)")
ax6c = Axis(fig[6, 3], xlabel = "Time (ms)", ylabel = "Gating Variable", title = "Na Inactivation (h)")
ax6d = Axis(fig[6, 4], xlabel = "Time (ms)", ylabel = "Gating Variable", title = "K Activation (n)")

# Plot all series
lines!(ax1a, tseries, atp_series, color = :blue, linewidth = 2)
lines!(ax1b, tseries, adp_series, color = :orange, linewidth = 2)
lines!(ax1c, tseries, amp_series, color = :green, linewidth = 2)
lines!(ax1d, tseries, ado_series, color = :red, linewidth = 2)
lines!(ax1e, tseries, pi_series, color = :purple, linewidth = 2)
lines!(ax1f, tseries, ribose1p_series, color = :magenta, linewidth = 2)

lines!(ax2a, tseries, pcr_series, color = :darkblue, linewidth = 2, label = "PCr")
lines!(ax2a, tseries, cr_series, color = :darkorange, linewidth = 2, label = "Cr")
lines!(ax2b, tseries, glucose_series, color = :brown, linewidth = 2)
lines!(ax2c, tseries, pep_series, color = :darkgreen, linewidth = 2)
lines!(ax2d, tseries, bpg_series, color = :darkred, linewidth = 2)
lines!(ax2e, tseries, pyruvate_series, color = :yellow, linewidth = 2)

lines!(ax3a, tseries, p3g_series, color = :brown, linewidth = 2)
lines!(ax3b, tseries, acetyl_coa_series, color = :olive, linewidth = 2)
lines!(ax3c, tseries, coa_series, color = :teal, linewidth = 2)

lines!(ax4a, tseries, nad_series, color = :forestgreen, linewidth = 2, label = "NAD")
lines!(ax4a, tseries, nadh_series, color = :hotpink, linewidth = 2, label = "NADH")
lines!(ax4b, tseries, fad_series, color = :khaki, linewidth = 2, label = "FAD")
lines!(ax4b, tseries, fadh2_series, color = :lavender, linewidth = 2, label = "FADH2")

lines!(ax5a, tseries, oxygen_series, color = :cyan, linewidth = 2, label = "O2")
lines!(ax5a, tseries, co2_series, color = :gray, linewidth = 2, label = "CO2")
lines!(ax5a, tseries, h2o_series, color = :lightblue, linewidth = 2, label = "H2O")

# Plot voltage and gating variables
lines!(ax6a, tseries, v_series, color = :black, linewidth = 3)
lines!(ax6b, tseries, m_series, color = :red, linewidth = 2)
lines!(ax6c, tseries, h_series, color = :blue, linewidth = 2)
lines!(ax6d, tseries, n_series, color = :green, linewidth = 2)

fig
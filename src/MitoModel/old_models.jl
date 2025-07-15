using DifferentialEquations
using GLMakie

#%% One-liner auxiliary functions

# Hodgkin-Huxley gating variables
alpha_m(Vm) = 0.1 * (Vm + 40.0) / (1.0 - exp(-(Vm + 40.0) / 10.0))
beta_m(Vm) = 4.0 * exp(-(Vm + 65.0) / 18.0)
alpha_h(Vm) = 0.07 * exp(-(Vm + 65.0) / 20.0)
beta_h(Vm) = 1.0 / (1.0 + exp(-(Vm + 35.0) / 10.0))
alpha_n(Vm) = 0.01 * (Vm + 55.0) / (1.0 - exp(-(Vm + 55.0) / 10.0))
beta_n(Vm) = 0.125 * exp(-(Vm + 65.0) / 80.0)

# Steady-state gating variables
m_infinity(Vm) = alpha_m(Vm) / (alpha_m(Vm) + beta_m(Vm))
h_infinity(Vm) = alpha_h(Vm) / (alpha_h(Vm) + beta_h(Vm))
n_infinity(Vm) = alpha_n(Vm) / (alpha_n(Vm) + beta_n(Vm))

# Applied current stimulation
applied_current(t, stim_start, stim_end, stim_amplitude) = (stim_start < t < stim_end) ? stim_amplitude : 0.0
pump_current(Vm, atp, I_pump_max, E_pump_target, K_m_ATP_pump) = I_pump_max * (atp / (atp + K_m_ATP_pump)) * max(0.0, (E_pump_target - Vm) / 30.0)
pump_atp_consumption(Vm, atp, I_pump_max, E_pump_target, K_m_ATP_pump, atp_per_current) = pump_current(Vm, atp, I_pump_max, E_pump_target, K_m_ATP_pump) * atp_per_current

# Metabolic functions
glucose_uptake_rate(glucose, k_glucose_uptake, K_m_glucose) = k_glucose_uptake * glucose / (K_m_glucose + glucose)
oxidative_phosphorylation_rate(adp, oxygen, V_max_ox_phos, K_m_oxygen) = V_max_ox_phos * adp * oxygen * adp / (adp + 0.01) * oxygen / (oxygen + 0.001) / ((K_m_oxygen + oxygen) * (1.0 + adp))
atp_synthesis_base_rate(adp, P, k_atp_synth) = k_atp_synth * adp * P * adp / (adp + 0.01) * P / (P + 0.01)
activity_factor(Vm) = 1.0 + 0.2 * max(0.0, (Vm + 70.0) / 40.0)
substrate_influx(current_level, max_level, diffusion_rate) = diffusion_rate * (max_level - current_level)
protection_factor(concentration, threshold = 0.001) = concentration / (concentration + threshold)

# Main ODE function
function mito_model_ode!(du, u, p, t)
    # Unpack variables
    atp, adp, amp, ado, P, Vm, glucose, oxygen = u
    
    # Unpack parameters
    k_atp_adp, k_adp_amp, k_amp_ado, k_adk, k_ak, k_atp_synth = p[1:6]
    g_LEAK, g_NaV, g_KV = p[7:9]
    E_LEAK, E_Na, E_K, C_m = p[10:13]
    stim_start, stim_end, stim_amplitude = p[14:16]
    
    # Na-K-ATPase pump parameters
    I_pump_max, E_pump_target, K_m_ATP_pump, atp_per_current = p[17:20]
    
    # Glucose and oxygen parameters
    glucose_max, K_m_glucose, k_glucose_uptake, V_max_ox_phos, K_m_oxygen = p[21:25]
    glucose_diffusion_rate, oxygen_max, oxygen_diffusion_rate = p[26:28]
    ATP_yield_glucose, ATP_yield_oxygen = p[29:30]

    # Gating variables using auxiliary functions
    m_steady = m_infinity(Vm)
    h_steady = h_infinity(Vm)
    n_steady = n_infinity(Vm)

    # Na-K-ATPase pump - using your function signatures
    I_pump = pump_current(Vm, atp, I_pump_max, E_pump_target, K_m_ATP_pump)
    pump_atp_consumed = pump_atp_consumption(Vm, atp, I_pump_max, E_pump_target, K_m_ATP_pump, atp_per_current)

    # Currents (including pump)
    I_LEAK = g_LEAK * (Vm - E_LEAK)
    I_NaV = g_NaV * m_steady^3 * h_steady * (Vm - E_Na)
    I_KV = g_KV * n_steady^4 * (Vm - E_K)
    # Applied current using auxiliary function
    I_APP = applied_current(t, stim_start, stim_end, stim_amplitude)

    # Metabolic rates using auxiliary functions
    glucose_uptake = glucose_uptake_rate(glucose, k_glucose_uptake, K_m_glucose)
    ox_phos_rate = oxidative_phosphorylation_rate(adp, oxygen, V_max_ox_phos, K_m_oxygen)
    base_atp_synth = atp_synthesis_base_rate(adp, P, k_atp_synth)
    
    # Enhanced ATP synthesis
    enhanced_atp_synth = base_atp_synth + 
                        ATP_yield_glucose * glucose_uptake + 
                        ATP_yield_oxygen * ox_phos_rate
    
    # Oxygen consumption
    act_factor = activity_factor(Vm)
    oxygen_consumption = 0.03 * act_factor + 0.15 * ox_phos_rate

    # ODEs with realistic pump ATP consumption
    du[1] = -k_atp_adp*atp - k_adk*atp*ado + k_ak*adp^2 + enhanced_atp_synth - pump_atp_consumed
    du[2] = k_atp_adp*atp - k_adp_amp*adp + k_adk*atp*ado - 2*k_ak*adp^2 - enhanced_atp_synth
    du[3] = k_adp_amp*adp - k_amp_ado*amp + k_adk*atp*ado + k_ak*adp^2
    du[4] = k_amp_ado*amp - k_adk*atp*ado
    du[5] = k_atp_adp*atp + k_adp_amp*adp - enhanced_atp_synth + pump_atp_consumed
    du[6] = -(I_LEAK + I_NaV + I_KV - I_pump - I_APP) / C_m
    
    # Substrate dynamics
    glucose_influx = substrate_influx(glucose, glucose_max, glucose_diffusion_rate)
    oxygen_influx = substrate_influx(oxygen, oxygen_max, oxygen_diffusion_rate)
    
    glucose_protection = protection_factor(glucose)
    oxygen_protection = protection_factor(oxygen)
    
    du[7] = glucose_influx - glucose_uptake * glucose_protection
    du[8] = oxygen_influx - oxygen_consumption * oxygen_protection

    return du
end

# Parameters with realistic Na-K-ATPase pump
p = [
    # Biochemical parameters
    0.00005,  # k_atp_adp (s⁻¹)
    0.00002,  # k_adp_amp (s⁻¹)
    0.002,    # k_amp_ado (s⁻¹)
    0.005,    # k_adk (mM⁻¹s⁻¹)
    0.010,    # k_ak (mM⁻¹s⁻¹)
    2.0,      # k_atp_synth (mM⁻¹s⁻¹)
    
    # Electrical parameters
    5e-5,     # g_LEAK (S/cm²)
    10.0,     # g_NaV (S/cm²)
    2.0,      # g_KV (S/cm²)
    -54.4,    # E_LEAK (mV)
    55.0,     # E_Na (mV)
    -77.0,    # E_K (mV)
    1.0,      # C_m (μF/cm²)
    
    # Stimulation parameters
    5.0,      # stim_start (ms)
    105.0,      # stim_end (ms)
    10.0,     # stim_amplitude (μA/cm²)
    
    # Na-K-ATPase pump parameters
    1.5,      # I_pump_max (μA/cm²) - maximum pump current
    -85.0,    # E_pump_target (mV) - target/equilibrium potential
    0.1,      # K_m_ATP_pump (mM) - ATP half-saturation for pump
    0.01,     # atp_per_current (mM*cm²/μA/s) - ATP consumed per unit current
    
    # Metabolic parameters
    10.0,     # glucose_max (mM)
    5.0,      # K_m_glucose (mM)
    1.2,      # k_glucose_uptake (mM/s)
    2.0,      # V_max_ox_phos (mM/s)
    0.05,     # K_m_oxygen (mM)
    0.20,     # glucose_diffusion_rate (s⁻¹)
    0.3,      # oxygen_max (mM)
    0.35,     # oxygen_diffusion_rate (s⁻¹)
    3.0,      # ATP_yield_glucose
    8.0       # ATP_yield_oxygen
]

# Initial conditions
u0 = [10.0, 0.0, 0.0, 0.0, 0.0, -65.0, 10.0, 0.3]  # [ATP, ADP, AMP, Ado, Pi, Vm, glucose, oxygen]
tspan = (0.0, 100.0)
prob = ODEProblem(mito_model_ode!, u0, tspan, p)
sol = solve(prob, Tsit5())

# Extract the data
tseries = LinRange(tspan[1], tspan[2], 1000)
atp_series = map(t -> sol(t)[1], tseries)
adp_series = map(t -> sol(t)[2], tseries)
amp_series = map(t -> sol(t)[3], tseries)
ado_series = map(t -> sol(t)[4], tseries)
pi_series = map(t -> sol(t)[5], tseries)
vm_series = map(t -> sol(t)[6], tseries)
glucose_series = map(t -> sol(t)[7], tseries)
oxygen_series = map(t -> sol(t)[8], tseries)

#%% Enhanced plotting with external legends
fig = Figure(size = (1600, 1000))

# Create layout
ga = fig[1, 1] = GridLayout()
gl = fig[1, 2] = GridLayout()

# Plot adenine nucleotides
ax1 = Axis(ga[1, 1], xlabel = "Time (ms)", ylabel = "Concentration (mM)", 
          title = "Adenine Nucleotides")

ax1a = Axis(ga[1, 1], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "ATP")
ax1b = Axis(ga[1, 2], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "ADP")
ax1c = Axis(ga[1, 3], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "AMP")
ax1d = Axis(ga[1, 4], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "Ado")
ax1e = Axis(ga[1, 5], xlabel = "Time (ms)", ylabel = "Concentration (mM)", title = "Pi")

lines!(ax1a, tseries, atp_series, color = :blue, linewidth = 2)
lines!(ax1b, tseries, adp_series, color = :orange, linewidth = 2)
lines!(ax1c, tseries, amp_series, color = :green, linewidth = 2)
lines!(ax1d, tseries, ado_series, color = :red, linewidth = 2)
lines!(ax1e, tseries, pi_series, color = :purple, linewidth = 2)

# Plot voltage
ax2 = Axis(ga[2, 1:4], xlabel = "Time (ms)", ylabel = "Voltage (mV)", 
          title = "Membrane Potential")
line_vm = lines!(ax2, tseries, vm_series, color = :black, linewidth = 2)

# Plot glucose dynamics
ax3 = Axis(ga[3, 1], xlabel = "Time (ms)", ylabel = "Concentration (mM)", 
          title = "Glucose Dynamics")
line_glucose = lines!(ax3, tseries, glucose_series, color = :brown, linewidth = 2)

# Plot oxygen dynamics
ax4 = Axis(ga[3, 3], xlabel = "Time (ms)", ylabel = "Concentration (mM)", 
          title = "Oxygen Dynamics")
line_oxygen = lines!(ax4, tseries, oxygen_series, color = :cyan, linewidth = 2)

# Plot ATP/ADP ratio
ax5 = Axis(ga[3, 2], xlabel = "Time (ms)", ylabel = "Ratio", 
          title = "ATP/ADP Ratio (Literature: 4-100)")
atp_adp_ratio = map(t -> sol(t)[1]/max(sol(t)[2], 1e-6), tseries)
line_ratio = lines!(ax5, tseries, atp_adp_ratio, color = :magenta, linewidth = 2)
hlines!(ax5, [4, 20, 50], color = :gray, linestyle = :dash, alpha = 0.5)
text!(ax5, 20, 45, text = "Typical range: 4-50", fontsize = 10)

# Plot energy charge
ax6 = Axis(ga[3, 4], xlabel = "Time (ms)", ylabel = "Energy Charge", 
          title = "Energy Charge ([ATP + 0.5×ADP]/[ATP + ADP + AMP])")
energy_charge = map(t -> (sol(t)[1] + 0.5*sol(t)[2]) / (sol(t)[1] + sol(t)[2] + sol(t)[3]), tseries)
line_energy = lines!(ax6, tseries, energy_charge, color = :darkgreen, linewidth = 2)
hlines!(ax6, [0.8, 0.9], color = :gray, linestyle = :dash, alpha = 0.5)
text!(ax6, 20, 0.85, text = "Healthy range: 0.8-0.9", fontsize = 10)

# External legends
Legend(gl[1, 1], lines1, labels1, "Adenine Nucleotides", framevisible = true)
Legend(gl[2, 1], [line_vm], ["Vm"], "Membrane Potential", framevisible = true)
Legend(gl[3, 1], [line_glucose], ["Glucose"], "Substrates", framevisible = true)
Legend(gl[4, 1], [line_oxygen], ["Oxygen"], "Oxygen", framevisible = true)
Legend(gl[5, 1], [line_ratio], ["ATP/ADP"], "Energy Ratios", framevisible = true)
Legend(gl[6, 1], [line_energy], ["Energy Charge"], "Energy Status", framevisible = true)

# Adjust column sizes
colsize!(fig.layout, 1, Relative(0.8))
colsize!(fig.layout, 2, Relative(0.2))

display(fig)

function pde_system!(du::Vector{T}, u::Vector{T}, p::Params2D, t; N = 1) where {T}
    if t % 10 == 0
        println("t = $t")
    end
    println("t = $t")

    # Split the solution vector into components
    ca = @view u[1:N]
    atp = @view u[N+1:2N]
    adp = @view u[2N+1:3N]
    amp = @view u[3N+1:4N]
    ado = @view u[4N+1:5N]
    P = @view u[5N+1:6N]  # Phosphate concentration
    Vm = @view u[6N+1:7N]  # Membrane voltage
    
    # Split the derivative vector
    dca = @view du[1:N]
    datp = @view du[N+1:2N]
    dadp = @view du[2N+1:3N]
    damp = @view du[3N+1:4N]
    dado = @view du[4N+1:5N]
    dP = @view du[5N+1:6N]  # Phosphate derivative
    dVm = @view du[6N+1:7N]  # Membrane voltage derivative
    
    # Convert masks to vectors for broadcasting
    cyto_mask_vec = p.cyto_mask
    mito_mask_vec = p.mito_mask  # Mitochondria mask is inverse of cytoplasm mask
    
    # Compute diffusion for all species (only in cytoplasm, with firm boundary at mitochondria)
    mul!(dca, p.L, ca)
    @. dca *= p.D_ca
    
    mul!(datp, p.L, atp)
    @. datp *= p.D_atp
    
    mul!(dadp, p.L, adp)
    @. dadp *= p.D_adp
    
    mul!(damp, p.L, amp)
    @. damp *= p.D_amp
    
    mul!(dado, p.L, ado)
    @. dado *= p.D_ado
    
    mul!(dP, p.L, P)
    @. dP *= p.D_p
    
    # Hodgkin-Huxley channel dynamics
    # ATP-dependent leak channel gating
    θ_ATP = hill_equation.(atp, p.Kd_ATP, p.hill_n)
    g_LEAK_eff = p.g_LEAK .+ p.gKATP_max .* θ_ATP
    
    # Leak current with ATP-dependent gating
    I_LEAK = g_LEAK_eff .* (Vm .- p.E_LEAK)
    
    # NaV channel (m³h)
    # Activation gate (m)
    α_m = 0.1 * (Vm .+ 40.0) ./ (1.0 .- exp.(-(Vm .+ 40.0) ./ 10.0))
    β_m = 4.0 * exp.(-(Vm .+ 65.0) ./ 18.0)
    m_inf = α_m ./ (α_m .+ β_m)
    τ_m = 1.0 ./ (α_m .+ β_m)
    
    # Inactivation gate (h)
    α_h = 0.07 * exp.(-(Vm .+ 65.0) ./ 20.0)
    β_h = 1.0 ./ (1.0 .+ exp.(-(Vm .+ 35.0) ./ 10.0))
    h_inf = α_h ./ (α_h .+ β_h)
    τ_h = 1.0 ./ (α_h .+ β_h)
    
    # NaV current
    I_NaV = p.g_NaV * m_inf.^3 .* h_inf .* (Vm .- p.E_Na)
    
    # KV channel (n⁴)
    # Activation gate (n)
    α_n = 0.01 * (Vm .+ 55.0) ./ (1.0 .- exp.(-(Vm .+ 55.0) ./ 10.0))
    β_n = 0.125 * exp.(-(Vm .+ 65.0) ./ 80.0)
    n_inf = α_n ./ (α_n .+ β_n)
    τ_n = 1.0 ./ (α_n .+ β_n)
    
    # KV current
    I_KV = p.g_KV * n_inf.^4 .* (Vm .- p.E_K)
    
    # Total current and voltage update
    I_APP = p.stim_start < t < p.stim_end ? p.stim_amplitude : 0.0
    I_total = I_LEAK .+ I_NaV .+ I_KV .- I_APP
    @. dVm = (-I_total / p.C_m) .* cyto_mask_vec
    
    # Calculate ATP consumption rate from ATP-dependent leak conductance
    leak_atp_rate = (p.gKATP_max .* θ_ATP) .* p.k_leak_atp .* cyto_mask_vec
    #println("leak_atp_rate = $(leak_atp_rate[1])")
    # Add cytoplasmic reaction terms
    # ATP → ADP conversion (in cytoplasm) - now includes leak current contribution
    @. datp -= (p.k_atp_adp * atp + leak_atp_rate) .* cyto_mask_vec
    @. dadp += (p.k_atp_adp * atp + leak_atp_rate) .* cyto_mask_vec
    @. dP += (p.k_atp_adp * atp + leak_atp_rate) .* cyto_mask_vec  # Release phosphate
    
    # ADP → AMP conversion (in cytoplasm)
    @. dadp -= p.k_adp_amp * adp .* cyto_mask_vec
    @. damp += p.k_adp_amp * adp .* cyto_mask_vec
    @. dP += p.k_adp_amp * adp .* cyto_mask_vec  # Release phosphate
    
    # AMP → Adenosine conversion (in cytoplasm)
    @. damp -= p.k_amp_ado * amp .* cyto_mask_vec
    @. dado += p.k_amp_ado * amp .* cyto_mask_vec
    
    # ADK reaction: Adenosine + ATP → AMP + ADP (in cytoplasm)
    @. datp -= p.k_adk * atp .* ado .* cyto_mask_vec  # ATP consumed
    @. dado -= p.k_adk * atp .* ado .* cyto_mask_vec  # Adenosine consumed
    @. damp += p.k_adk * atp .* ado .* cyto_mask_vec  # AMP produced
    @. dadp += p.k_adk * atp .* ado .* cyto_mask_vec  # ADP produced
    
    # Add mitochondrial reaction terms
    # Mitochondrial adenylate kinase reaction: 2ADP → ATP + AMP (only in mitochondria)
    @. datp += p.k_ak * mito_mask_vec .* adp.^2  # ATP produced
    @. dadp -= 2.0 * p.k_ak * mito_mask_vec .* adp.^2  # 2 ADP consumed
    @. damp += p.k_ak * mito_mask_vec .* adp.^2  # AMP produced
    
    # ATP synthesis from ADP + P (only in mitochondria)
    @. datp += p.k_atp_synth * mito_mask_vec .* adp .* P  # ATP produced
    @. dadp -= p.k_atp_synth * mito_mask_vec .* adp .* P  # ADP consumed
    @. dP -= p.k_atp_synth * mito_mask_vec .* adp .* P    # Phosphate consumed
    
    return nothing
end


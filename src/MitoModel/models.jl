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
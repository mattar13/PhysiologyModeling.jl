function pde_system!(du::Vector{Float64}, u::Vector{Float64}, p::Params2D, t; N = 1)
    if t % 10 == 0
        println("t = $t")
    end
    
    # Split the solution vector into components
    ca = @view u[1:N]
    atp = @view u[N+1:2N]
    adp = @view u[2N+1:3N]
    amp = @view u[3N+1:4N]
    ado = @view u[4N+1:5N]
    P = @view u[5N+1:6N]  # Phosphate concentration
    
    # Split the derivative vector
    dca = @view du[1:N]
    datp = @view du[N+1:2N]
    dadp = @view du[2N+1:3N]
    damp = @view du[3N+1:4N]
    dado = @view du[4N+1:5N]
    dP = @view du[5N+1:6N]  # Phosphate derivative
    
    # Convert masks to vectors for broadcasting
    cyto_mask_vec = p.cyto_mask
    mito_mask_vec = p.mito_mask  # Mitochondria mask is inverse of cytoplasm mask
    
    # Initialize derivatives to zero
    # fill!(du, 0.0)
    
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
    
    # Add cytoplasmic reaction terms
    # ATP → ADP conversion (in cytoplasm)
    @. datp -= p.k_atp_adp * atp .* cyto_mask_vec
    @. dadp += p.k_atp_adp * atp .* cyto_mask_vec
    @. dP += p.k_atp_adp * atp .* cyto_mask_vec  # Release phosphate
    
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
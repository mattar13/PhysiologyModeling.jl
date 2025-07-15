# Create metabolic reaction network
metabolic_network = @reaction_network begin
    # ATP hydrolysis to ADP + Pi
    k_atp_adp, ATP --> ADP + Pi
    
    # ADP conversion to AMP + Pi  
    k_adp_amp, ADP --> AMP + Pi
    
    # AMP conversion to Adenosine + Pi
    k_amp_ado, AMP --> Ado + Pi
    
    # Adenylate kinase reaction: ATP + AMP <--> 2 ADP
    k_adk_forward, ATP + AMP --> 2*ADP
    k_adk_reverse, 2*ADP --> ATP + AMP
    
    # ATP synthesis from ADP + Pi (mitochondrial)
    k_atp_synth, ADP + Pi --> ATP
    
    # Adenosine salvage pathway
    k_ado_kinase, Ado + ATP --> AMP + ADP
    k_ado_salvage, Ado + R1P --> AMP
    
    # Glycolysis: Glucose → 2 Pyruvate + 2 ATP
    k_glycolysis, GLU + 2*ADP + 2*Pi --> 2*Pyruvate + 2*ATP
    
    # Pyruvate oxidation: Pyruvate → Acetyl-CoA  
    k_pyruvate_dehydrogenase, Pyruvate + CoA --> Acetyl_CoA + CO2
    
    # TCA cycle: Acetyl-CoA → 3 NADH + 1 FADH2 + 1 ATP
    k_tca_cycle, Acetyl_CoA + 3*NAD + FAD + ADP + Pi --> 3*NADH + FADH2 + ATP + 2*CO2
    
    # Electron transport: NADH/FADH2 → ATP
    k_complex_I, 2*NADH + 6*ADP + 6*Pi + O2 --> 2*NAD + 6*ATP + 2*H2O
    k_complex_II, 2*FADH2 + 4*ADP + 4*Pi + O2 --> 2*FAD + 4*ATP + 2*H2O
    
    # Formation from glucose (simplified glycolysis)
    k_glucose_to_pep, GLU + Pi --> PEP
    k_glucose_to_bpg, GLU + Pi --> BPG
    
    # Bootstrap reactions - convert AMP to ADP using high-energy intermediates
    k_pyruvate_kinase, AMP + PEP --> ADP + Pyruvate
    k_pgk, AMP + BPG --> ADP + P3G
    
    # Glucose transport
    kGLUT(GLU, k_glut_max), 0 --> GLU

    # Oxygen replenishment
    k_oxygen_supply, 0 --> O2
        
    # Basal consumption
    k_glucose_basal, GLU --> 0
    k_oxygen_basal, O2 --> 0
    
    # Creatine kinase system
    k_creatine_kinase_forward, ATP + Cr --> ADP + PCr
    k_creatine_kinase_reverse, PCr + ADP --> ATP + Cr
    
    # Na/K ATPase consumption (ATP-dependent)
    I_NaK_ATPase(ATP, I_NaK_max, K_ATP_NaK), ATP --> ADP + Pi
end

# Convert to ODESystem and get metabolic equations
metabolic_sys = convert(ODESystem, metabolic_network)
metabolic_eqs = equations(metabolic_sys)

# Create HH equations manually
hh_eqs = [
    D(m) ~ (m_inf(V) - m) / τ_m(V),
    D(h) ~ (h_inf(V) - h) / τ_h(V),
    D(n) ~ (n_inf(V) - n) / τ_n(V),
    D(V) ~ -(I_Na(V, m, h, g_Na_ATP(ATP, g_Na_max, K_ATP_Na), E_Na) + 
             I_K(V, n, g_K_ATP(ATP, g_K_max, K_ATP_K), E_K) + 
             I_leak(V, g_leak, E_leak) - 
             #I_NaK_ATPase(ATP, I_NaK_max, K_ATP_NaK) + 
             I_APP(t, stim_start, stim_end, I_amplitude)) / C_m
]

# Combine all equations
all_eqs = vcat(metabolic_eqs, hh_eqs)

# Create the complete system
@named complete_system = ODESystem(all_eqs, t)

# Complete the system
complete_system = structural_simplify(complete_system)

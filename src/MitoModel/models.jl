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
    
    # Adenosine salvage pathway
    k_ado_kinase, Ado + ATP --> AMP + ADP
    k_ado_salvage, Ado + Pi --> AMP
    
    # Glycolysis with feedback inhibition: Glucose → 2 Pyruvate + 2 ATP 
    # Rate is modulated by pH, pyruvate, ATP, alanine, and H+ levels
    k_glycolysis * glycolysis_inhibition(Lactate, Pyruvate, ATP, Alanine, H, K_pH, K_pyruvate, K_ATP_glyc, K_alanine, K_H), GLU + 2*ADP + 2*Pi + 2*NAD --> 2*Pyruvate + 2*ATP + 2*NADH

    # Pyruvate oxidation: Pyruvate → Acetyl-CoA  
    k_pyruvate_dehydrogenase, Pyruvate + CoA --> Acetyl_CoA + CO2
    
    # TCA cycle: Acetyl-CoA → 3 NADH + 1 FADH2 + 1 ATP + 2 CO2 + CoA
    k_tca_cycle, Acetyl_CoA + 3*NAD + FAD + ADP + Pi --> 3*NADH + FADH2 + ATP + 2*CO2 + CoA

    # Electron transport: NADH → 2.5 ATP (but using integers: 2 NADH → 5 ATP)
    k_complex_I, 2*NADH + 5*ADP + 5*Pi + O2 --> 2*NAD + 5*ATP + H2O

    # Electron transport: FADH2 → 1.5 ATP (but using integers: 2 FADH2 → 3 ATP)  
    k_complex_II, 2*FADH2 + 3*ADP + 3*Pi + O2 --> 2*FAD + 3*ATP + H2O
    
    # === Anaerobic Metabolism ===
    k_lactate_production, Pyruvate + NADH --> Lactate + NAD
    k_LDH, Lactate + NAD --> Pyruvate + NADH
    
    # === Alanine Metabolism ===
    # Bidirectional alanine-pyruvate conversion
    k_alanine_synthesis, Pyruvate --> Alanine        # High pyruvate → alanine
    k_alanine_breakdown, Alanine --> Pyruvate        # Can reverse when pyruvate drops
    
    # Alanine disposal (export or utilization)
    k_alanine_disposal, Alanine --> 0               # Represents export or protein synthesis
    k_lactate_disposal, Lactate --> 0
    
    # Glucose transport
    kGLUT(GLU, k_glut_supply, k_glut_extracellular), 0 --> GLU

    # Oxygen replenishment
    kOXYGEN(O2, k_oxygen_supply, k_oxygen_extracellular), 0 --> O2
        
    # Basal consumption
    k_glut_basal, GLU --> 0
    k_oxygen_basal, O2 --> 0
    
    # Creatine kinase system
    k_creatine_kinase_forward, ATP + Cr --> ADP + PCr
    k_creatine_kinase_reverse, PCr + ADP --> ATP + Cr
    
    # === CARBONIC ACID EQUILIBRIUM ===
    k_hydration,    CO2 + H2O       --> H2CO3
    k_dehydration,  H2CO3           --> CO2 + H2O
    k_dissociation, H2CO3           --> HCO3 + 2*H
    k_assoc,        HCO3 + 2*H        --> H2CO3
    
    #CO2 and H2O disposal
    k_co2_disposal, CO2 --> 0
    k_h2o_disposal, H2O --> 0

    # === ION FLUXES AND PUMP ===
    k_pump_max*a*d, ATP --> ADP + Pi    # 1 ATP per cycle, now product of activation and drive
        
    # Add the voltage equations for the HH model here
    @parameters C_m g_Na g_K g_leak E_Na E_K E_leak stim_start stim_end I_amplitude I_pump_max a_ATP_half a_V_half a_V_slope
    @equations begin
        # Remove dynamic reversal potentials and any equations involving Na_i, K_i, Na_o, K_o
        # HH equations with dynamic E_Na, E_K
        D(V) ~ (
            - g_Na * m^3 * h * (V - a*E_Na)
            - g_K * n^4 * (V - a*E_K)
            - g_leak * (V - E_leak)
            #+ I_pump_max * a * d
            + I_amplitude * min(1.0, max(0.0, (t - stim_start))) * min(1.0, max(0.0, (stim_end - t)))  # Step function from stim_start to stim_end
        ) / C_m
        
        # Direct HH equations without function calls
        D(h) ~ (0.07 * exp(-(V + 65) / 20) * (1 - h) - (1 / (1 + exp(-(V + 35) / 10))) * h)
        D(m) ~ (0.1 * (V + 40) / (1 - exp(-(V + 40) / 10)) * (1 - m) - 4 * exp(-(V + 65) / 18) * m)
        D(n) ~ (0.01 * (V + 55) / (1 - exp(-(V + 55) / 10)) * (1 - n) - 0.125 * exp(-(V + 65) / 80) * n)
        # Pump activation: now split into metabolic (a) and voltage (d) components
        D(a) ~ (ATP/(ATP + a_ATP_half)) - a
        D(d) ~ (1/(1 + exp(-(V - a_V_half)/a_V_slope))) - d
        
    end
end
@variables t
@species ATP(t) ADP(t) AMP(t) Ado(t) Pi(t) GLU(t) O2(t) PCr(t) Cr(t) CO2(t) H2O(t) PEP(t) BPG(t) Pyruvate(t) P3G(t) Acetyl_CoA(t) CoA(t) NAD(t) NADH(t) FAD(t) FADH2(t) Lactate(t) Alanine(t)
# Add ionic species
@species Na_i(t) K_i(t)
@variables V(t) m(t) h(t) n(t)
@variables a(t) d(t)

# Metabolic parameters
@parameters k_atp_adp k_adp_amp k_amp_ado k_adk_forward k_adk_reverse k_ado_kinase k_ado_salvage
@parameters k_glut_supply k_glut_extracellular 
@parameters k_oxygen_supply k_oxygen_extracellular
@parameters k_glut_basal k_oxygen_basal

@parameters k_glucose_to_pep k_glucose_to_bpg k_pyruvate_kinase k_pgk k_p3g_utilization
@parameters k_glycolysis k_pyruvate_dehydrogenase k_tca_cycle k_complex_I k_complex_II
@parameters k_creatine_kinase_forward k_creatine_kinase_reverse
@parameters k_lactate_production k_LDH k_lactate_disposal
@parameters K_pH K_pyruvate K_ATP_glyc k_alanine_synthesis K_alanine
@parameters k_alanine_breakdown k_alanine_disposal
# Hodgkin-Huxley parameters
@parameters C_m g_Na g_K g_leak E_Na E_K E_leak K_ATP_NaK I_NaK_max stim_start stim_end I_amplitude
# Pump activation parameters
@parameters a_V_half a_V_slope a_ATP_half I_pump_max k_pump_max

D = Differential(t)

# Parameters for the combined system
params = (
    # ==== Metabolic parameters ====
    k_atp_adp => 0.02,           # | 0.1 | ATP -> ADP + Pi | ATP hydrolysis rate
    k_adp_amp => 0.0005,        # | 0.0005 | ADP -> AMP + Pi | ADP to AMP conversion
    k_amp_ado => 0.002,         # | 0.002 | AMP -> Ado + Pi | AMP to adenosine conversion
    k_adk_forward => 0.01,      # | 0.01 | ATP + AMP -> 2 ADP | Adenylate kinase (forward)
    k_adk_reverse => 0.02,      # | 0.02 | 2 ADP -> ATP + AMP | Adenylate kinase (reverse)
    k_ado_kinase => 0.005,      # | 0.005 | Ado + ATP -> AMP + ADP | Adenosine kinase
    k_ado_salvage => 0.001,     # | 0.001 | Ado + Pi -> AMP | Adenosine salvage pathway

    # ==== Glucose and oxygen supply ====
    k_glut_supply => 2.0,           # | 2.0 | Glucose supply rate | GLU
    k_glut_extracellular => 10.0,   # | 10.0 | Extracellular glucose | GLU
    k_oxygen_supply => 2.0,         # | 1.0 | Oxygen supply rate | O2
    k_oxygen_extracellular => 3.0,  # | 0.3 | Extracellular oxygen | O2
    
    # ==== Basal consumption ====
    k_glut_basal => 0.0,        # | 0.0 | Basal glucose consumption | GLU
    k_oxygen_basal => 0.0,      # | 0.0 | Basal oxygen consumption | O2

    # ==== Glycolysis ==== #
    k_glycolysis => 1.0,                # | 2.0 | Glycolysis rate | GLU, Pyruvate, ATP, NAD, NADH
    k_pyruvate_dehydrogenase => 0.3,    # | 0.3 | Pyruvate to Acetyl-CoA | Pyruvate, Acetyl_CoA
    
    # ==== TCA cycle ==== #
    k_tca_cycle => 0.05,                 # | 0.3 | TCA cycle rate | Acetyl_CoA, NADH, ATP
    k_complex_I => 2.0,                 # | 5.0 | ETC complex I rate | NADH, ATP, O2
    k_complex_II => 1.0,                # | 3.0 | ETC complex II rate | FADH2, ATP, O2
    
    # ==== Glucose to PEP and BPG ==== #
    k_glucose_to_pep => 0.2,            # | 0.2 | Glucose to PEP | GLU, PEP
    k_glucose_to_bpg => 0.15,           # | 0.15 | Glucose to BPG | GLU, BPG
    
    # ==== Pyruvate kinase system ==== #
    k_pyruvate_kinase => 0.8,           # | 0.8 | AMP + PEP -> ADP + Pyruvate | AMP, PEP, ADP, Pyruvate
    k_pgk => 0.6,                       # | 0.6 | AMP + BPG -> ADP + P3G | AMP, BPG, ADP, P3G
    k_p3g_utilization => 0.5,           # | 0.5 | P3G -> PEP | P3G utilization back to PEP

    # ==== Creatine kinase system ==== #
    k_creatine_kinase_forward => 1.0, # | 5.0 | ATP + Cr -> ADP + PCr | Creatine kinase (forward)
    k_creatine_kinase_reverse => 5.0, # | 2.0 | PCr + ADP -> ATP + Cr | Creatine kinase (reverse)
    
    # ==== Anaerobic Metabolism ====
    k_lactate_production => 0.5,   # | 0.5 | Anaerobic lactate production | Pyruvate, NADH, Lactate, NAD
    k_LDH => 1.0,                  # | 1.0 | Lactate dehydrogenase | Lactate, NADH, Pyruvate, NAD
    k_lactate_disposal => 0.1,     # | 0.5 | Lactate disposal | Lactate
    # ==== Feedback Inhibition Parameters ====
    K_pH => 2.0,                # | 2.0 | Lactate concentration for half-maximal pH inhibition | Lactate
    K_pyruvate => 1.0,           # | 1.0 | Pyruvate concentration for half-maximal inhibition | Pyruvate
    K_ATP_glyc => 0.5,           # | 0.5 | ATP concentration for glycolysis inhibition | ATP
    k_alanine_synthesis => 0.1,  # | 0.1 | Pyruvate to alanine conversion rate | Pyruvate, Alanine
    K_alanine => 0.5,            # | 0.5 | Alanine inhibition constant | Alanine
    
    # ==== Alanine Metabolism ====
    k_alanine_breakdown => 0.05, # | 0.05 | Alanine to pyruvate conversion rate | Alanine, Pyruvate
    k_alanine_disposal => 0.02,  # | 0.02 | Alanine disposal/export rate | Alanine
    
    # ==== Hodgkin-Huxley parameters ====
    C_m => 1.0,                 # | 1.0 | Membrane capacitance | V
    g_Na => 120.0,              # | 120.0 | Na conductance | V, m, h
    g_K => 36.0,                # | 36.0 | K conductance | V, n
    g_leak => 0.3,              # | 0.3 | Leak conductance | V
    E_Na => 50.0,               # | 50.0 | Na reversal potential | V
    E_K => -77.0,               # | -77.0 | K reversal potential | V
    E_leak => -54.4,            # | -54.4 | Leak reversal potential | V

    # ==== Stimulus parameters ====
    stim_start => 500.0,        # | 100.0 | Stimulus start time (ms) | V
    stim_end => 2500.0,          # | 500.0 | Stimulus end time (ms) | V
    I_amplitude => 10.0,        # | 10.0 | Stimulus amplitude | V

    # ==== Pump activation parameters ====
    a_V_half => -30.0,          # | -30.0 | Pump V half-activation | a, V
    a_V_slope => 10.0,          # | 10.0 | Pump V slope | a, V
    a_ATP_half => 0.3,          # | 0.5 | Pump ATP half-activation | a, ATP
    I_pump_max => 2.0,          # | 2.0 | Max pump current | V, a
    k_pump_max => 0.1           # | 1.0 | Max pump rate | ATP, ADP, Pi, a
)

# Initial conditions
u0 = [
    # === ADENINE NUCLEOTIDES === (sources: NCBI, PNAS studies)
    ATP => 5.0,      # 1-10 mM typical range, brain uses ~25% of body's ATP
    ADP => 0.5,    # Usually 10-20% of ATP concentration  
    AMP => 0.1,      # Very low under normal conditions (~2-5% of ATP)
    
    # === PHOSPHATE SYSTEM === (sources: Nature Communications, PMC studies)
    Pi => 50.0,       # 1-5 mM typical intracellular, CSF has ~100-fold lower
    PCr => 15.0,     # ~10-20 mM in brain tissue, major energy buffer
    Cr => 5.0,       # ~5-10 mM, total Cr+PCr ~20-25 mM in brain
    
    # === SUBSTRATES === (sources: Physiological studies)
    GLU => 0.0, #5.0,      # 2-10 mM typical brain glucose
    O2 => 0.0, #0.2,       # Low resting, rapidly consumed
    
    # === PURINE SALVAGE === (sources: Biochemical literature)
    Ado => 0.1,      # Micromolar range, neuroprotective
    
    # === METABOLIC INTERMEDIATES === (sources: Glycolysis studies)
    PEP => 0.0,      # Rapidly consumed in glycolysis
    BPG => 0.0,      # Rapidly consumed in glycolysis
    Pyruvate => 0.0, # Product formation
    P3G => 0.0,      # Product formation
    
    # === OXIDATIVE METABOLISM === (sources: Mitochondrial studies)
    Acetyl_CoA => 0.0,  # Rapidly consumed in TCA cycle
    CoA => 2.0,         # ~1-5 mM coenzyme A pool
    NAD => 10.0,        # ~5-15 mM total NAD pool
    NADH => 0.5,        # ~5-10% of total NAD under resting conditions
    FAD => 2.0,         # ~1-5 mM flavin pool
    FADH2 => 0.1,       # Small fraction of total FAD
    Lactate => 0.0,     # Anaerobic end product, starts at zero
    Alanine => 0.0,     # Amino acid product, starts at zero
    
    # === WASTE PRODUCTS ===
    CO2 => 0.0,      # Rapidly cleared
    H2O => 0.0,      # Abundant, not limiting
    
    # === ELECTRICAL VARIABLES === (validated HH values)
    V => -65.0,      # Standard resting potential
    h => 0.646593952660412,     # Na+ inactivation at rest (calculated)  
    m => 0.04446684116772448,     # Na+ activation at rest (calculated)
    n => 0.2953812119438617,      # K+ activation at rest (calculated)
    a => 1.0,
    d => 1.0
]
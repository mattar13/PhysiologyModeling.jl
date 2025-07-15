@variables t
@species ATP(t) ADP(t) AMP(t) Ado(t) Pi(t) GLU(t) O2(t) R1P(t) PCr(t) Cr(t) CO2(t) H2O(t) PEP(t) BPG(t) Pyruvate(t) P3G(t) Acetyl_CoA(t) CoA(t) NAD(t) NADH(t) FAD(t) FADH2(t)
# Add ionic species
@species Na_i(t) K_i(t)
@variables V(t) m(t) h(t) n(t)
@variables a(t)

# Metabolic parameters
@parameters k_atp_adp k_adp_amp k_amp_ado k_adk_forward k_adk_reverse k_ado_kinase k_ado_salvage
@parameters k_glut_supply k_glut_extracellular 
@parameters k_oxygen_supply k_oxygen_extracellular
@parameters k_glut_basal k_oxygen_basal

@parameters k_glucose_to_pep k_glucose_to_bpg k_pyruvate_kinase k_pgk
@parameters k_glycolysis k_pyruvate_dehydrogenase k_tca_cycle k_complex_I k_complex_II
@parameters k_creatine_kinase_forward k_creatine_kinase_reverse
# Hodgkin-Huxley parameters
@parameters C_m g_Na g_K g_leak E_Na E_K E_leak K_ATP_NaK I_NaK_max stim_start stim_end I_amplitude
# Pump activation parameters
@parameters a_V_half a_V_slope a_ATP_half I_pump_max k_pump_max

D = Differential(t)

# Parameters for the combined system
params = (
    # Metabolic parameters
    k_atp_adp => 0.001,
    k_adp_amp => 0.0005,
    k_amp_ado => 0.002,
    k_adk_forward => 0.01,
    k_adk_reverse => 0.02,
    k_ado_kinase => 0.005,
    k_ado_salvage => 0.001,

    #How glucose and oxygen are supplied to the cell
    k_glut_supply => 0.0,
    k_glut_extracellular => 10.0,
    k_oxygen_supply => 0.0,
    k_oxygen_extracellular => 0.3,
    
    #Basal consumption
    k_glut_basal => 0.0,
    k_oxygen_basal => 0.0,

    #Glycolysis
    k_glycolysis => 0.5,
    k_pyruvate_dehydrogenase => 0.3,
    k_tca_cycle => 0.3,
    k_complex_I => 2.5,
    k_complex_II => 1.5,
    k_glut_basal => 0.0,
    k_oxygen_basal => 0.0,
    k_glucose_to_pep => 0.2,
    k_glucose_to_bpg => 0.15,
    k_pyruvate_kinase => 0.8,
    k_pgk => 0.6,

    #Creatine kinase system
    k_creatine_kinase_forward => 0.3,
    k_creatine_kinase_reverse => 0.1,
    
    # Hodgkin-Huxley parameters
    C_m => 1.0,
    g_Na => 120.0,
    g_K => 36.0,
    g_leak => 0.3,
    E_Na => 50.0,
    E_K => -77.0,
    E_leak => -54.4,

    #Stimulus parameters    
    stim_start => 5.0,
    stim_end => 100.0,
    I_amplitude => 20.0,

    # Pump activation parameters
    a_V_half => -30.0,
    a_V_slope => 10.0,
    a_ATP_half => 0.5,
    I_pump_max => 0.5,
    k_pump_max => 10.0
)

# Initial conditions
u0 = [
    # === ADENINE NUCLEOTIDES === (sources: NCBI, PNAS studies)
    ATP => 3.0,      # 1-10 mM typical range, brain uses ~25% of body's ATP
    ADP => 0.5,      # Usually 10-20% of ATP concentration  
    AMP => 0.1,      # Very low under normal conditions (~2-5% of ATP)
    
    # === PHOSPHATE SYSTEM === (sources: Nature Communications, PMC studies)
    Pi => 2.0,       # 1-5 mM typical intracellular, CSF has ~100-fold lower
    PCr => 0.0,     # ~10-20 mM in brain tissue, major energy buffer
    Cr => 15.0,       # ~5-10 mM, total Cr+PCr ~20-25 mM in brain
    
    # === SUBSTRATES === (sources: Physiological studies)
    GLU => 0.0, #5.0,      # 2-10 mM typical brain glucose
    O2 => 0.0, #0.2,       # Low resting, rapidly consumed
    
    # === PURINE SALVAGE === (sources: Biochemical literature)
    Ado => 0.1,      # Micromolar range, neuroprotective
    R1P => 0.5,      # Ribose-1-phosphate for salvage pathway
    
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
    
    # === WASTE PRODUCTS ===
    CO2 => 0.0,      # Rapidly cleared
    H2O => 0.0,      # Abundant, not limiting
    
    # === ELECTRICAL VARIABLES === (validated HH values)
    V => -65.0,      # Standard resting potential
    h => 0.646593952660412,     # Na+ inactivation at rest (calculated)  
    m => 0.04446684116772448,     # Na+ activation at rest (calculated)
    n => 0.2953812119438617,      # K+ activation at rest (calculated)
    a => 1.0
]
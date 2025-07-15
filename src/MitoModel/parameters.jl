@variables t
@species ATP(t) ADP(t) AMP(t) Ado(t) Pi(t) GLU(t) O2(t) R1P(t) PCr(t) Cr(t) CO2(t) H2O(t) PEP(t) BPG(t) Pyruvate(t) P3G(t) Acetyl_CoA(t) CoA(t) NAD(t) NADH(t) FAD(t) FADH2(t)
@variables V(t) m(t) h(t) n(t)

# Metabolic parameters
@parameters k_atp_adp k_adp_amp k_amp_ado k_adk_forward k_adk_reverse k_atp_synth k_ado_kinase k_ado_salvage
@parameters k_glut_max k_oxygen_supply k_glycolysis k_pyruvate_dehydrogenase k_tca_cycle k_complex_I k_complex_II
@parameters k_glucose_basal k_oxygen_basal k_glucose_to_pep k_glucose_to_bpg k_pyruvate_kinase k_pgk
@parameters k_creatine_kinase_forward k_creatine_kinase_reverse

# Hodgkin-Huxley parameters
@parameters C_m g_Na_max g_K_max g_leak E_Na E_K E_leak K_ATP_Na K_ATP_K K_ATP_NaK I_NaK_max I_stim

D = Differential(t)

# Parameters for the combined system
params = (
    # Metabolic parameters
    k_atp_adp => 0.001,
    k_adp_amp => 0.0005,
    k_amp_ado => 0.002,
    k_adk_forward => 0.01,
    k_adk_reverse => 0.02,
    k_atp_synth => 0.5,
    k_ado_kinase => 0.005,
    k_ado_salvage => 0.001,
    k_glut_max => 0.0,
    k_oxygen_supply => 0.0,
    k_glycolysis => 0.5,
    k_pyruvate_dehydrogenase => 0.3,
    k_tca_cycle => 0.3,
    k_complex_I => 2.5,
    k_complex_II => 1.5,
    k_glucose_basal => 0.0,
    k_oxygen_basal => 0.0,
    k_glucose_to_pep => 0.2,
    k_glucose_to_bpg => 0.15,
    k_pyruvate_kinase => 0.8,
    k_pgk => 0.6,
    k_creatine_kinase_forward => 0.3,
    k_creatine_kinase_reverse => 0.1,
    
    # Hodgkin-Huxley parameters
    C_m => 1.0,
    g_Na_max => 120.0,
    g_K_max => 36.0,
    g_leak => 0.3,
    E_Na => 50.0,
    E_K => -77.0,
    E_leak => -54.387,
    K_ATP_Na => 2.0,
    K_ATP_K => 1.5,
    K_ATP_NaK => 0.5,
    I_NaK_max => 2.0,
    I_stim => 10.0
)

# Initial conditions
u0 = [
    Pi => 100.0,
    ATP => 0.0,
    ADP => 0.0,
    AMP => 0.0,
    Ado => 5.0,
    R1P => 1.0,
    GLU => 10.0,
    O2 => 0.0,
    PCr => 5.0,
    Cr => 0.0,
    CO2 => 0.0,
    H2O => 0.0,
    PEP => 0.0,
    BPG => 0.0,
    Pyruvate => 0.0,
    P3G => 0.0,
    Acetyl_CoA => 0.0,
    CoA => 10.0,
    NAD => 50.0,
    NADH => 0.0,
    FAD => 10.0,
    FADH2 => 0.0,
    V => -65.0,
    m => 0.0529,  # Initial m value at -65 mV
    h => 0.596,   # Initial h value at -65 mV
    n => 0.317    # Initial n value at -65 mV
]
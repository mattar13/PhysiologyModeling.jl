PC_p0_dict = Dict(
    "Cm" => 0.02,  # nF

    # Photocurrent parameters
    "a1" => 50,  # s^-1
    "a2" => 0.0003,  # s^-1
    "a3" => 0.03,  # s^-1
    "e" => 0.5,  # s^-1 μM^-1
    "b1" => 2.5,  # s^-1
    "s1" => 0.2,  # s^-1 μM^-1
    "s2" => 5,  # s^-1
    "J_max" => 5040,  # pA
    "Ttot" => 1000,  # μM
    "PDEtot" => 100,
    "cCa" => 50,  # s^-1
    "C0" => 0.1,  # μM
    "b" => 0.25,  # μM s^-1 pA^-1
    "k1" => 0.2,  # s^-1 μM^-1
    "k2" => 0.8,  # s^-1
    "eT" => 500,  # μM
    "V" => 0.4,  # s^-1
    "Kc" => 0.1,
    "Amax" => 65.6,  # μM s^-1
    "r" => 1.0,  # s^-1 μM^-1

    # Calcium-activated potassium current (IK(Ca)) parameters
    "gKCa" => 5.0,  # nS
    "EK" => -74,  # mV

    # Hyperpolarization-activated current (Ih) parameters
    "gh" => 3.0,  # nS
    "Eh" => -32,  # mV

    # Leakage current (IL) parameters
    "gL" => 0.35,  # nS
    "EL" => -77,  # mV

    # Delayed rectifier current (IKv) parameters
    "gKv" => 2.0,  # nS

    # Calcium current (ICa) parameters
    "gCa" => 0.7,  # nS
    "Ca_o" => 1600,  # μM

    # Calcium-activated chloride current (ICl(Ca)) parameters
    "gCl" => 2.0,  # nS
    "ECl" => -20  # mV

    # Intracellular calcium system parameters
    # Add any additional parameters here as needed
)

PC_u0_dict = Dict(
     "V" => -36.186,  # Membrane potential in mV
 
     # Photocurrent components
     "Rh" => 0,
     "Rhi" => 0,
     "Tr" => 0,
     "PDE" => 0,
     "cGMP" => 2.0,
 
     # Calcium and its buffer components
     "Ca" => 0.3,   # μM
     "Cab" => 34.88,  # μM
 
     # Ionic currents and their gating variables
     "mKCa" => 0.642,
     "C1" => 0.646,  # For Ih
     "C2" => 0.298,  # For Ih
     "O1" => 0.0517, # For Ih
     "O2" => 0.00398, # For Ih
     "O3" => 0.000115, # For Ih
     "mKv" => 0.430,
     "hKv" => 0.999,
     "mCa" => 0.436
)
 
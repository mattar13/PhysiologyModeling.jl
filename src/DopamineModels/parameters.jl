p = (
    Cm        = 0.5,      # pF
    gL        = 5.0,      # nS
    EL        = -70.0,    # mV
    gCa       = 1.5,      # nS
    ECa       = +120.0,   # mV
    Vhalf     = -25.0,    # mV
    kV        = 6.0,      # mV slope

    τCa       = 30.0,     # ms
    αCa       = 0.02,     # µM / (nA·ms)

    krel      = 0.05,     # µM / (µM·ms)
    kclear    = 0.001,     # ms⁻¹

    kon       = 0.02,     # µM⁻¹·ms⁻¹
    koff      = 0.01,     # ms⁻¹

    kG        = 0.04,     # µM·ms⁻¹
    kGdeg     = 0.005,    # ms⁻¹
    G50       = 0.2,      # µM
    nGi       = 2.0,      # Hill coeff.

    kACbasal  = 0.3,      # µM·ms⁻¹
    kcAMPdeg  = 0.005,    # ms⁻¹

    kPKA      = 0.03,     # µM·ms⁻¹
    kPKAdeg   = 0.01      # ms⁻¹
)

conds = (

)

settings = (
    tspan = (0.0, 500.0), # ms
    dt    = 0.5,          # ms
    nx    = 50,
    ny    = 50,
    dx    = 1.0,
    dy    = 1.0,
    D     = 0.3         # µM·ms⁻¹
)
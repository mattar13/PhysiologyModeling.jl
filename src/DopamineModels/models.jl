# ---------- ODE system ----------
function dopaminergic_autoreceptor!(du, u, p, t)
    (Cm, gL, EL, gCa, ECa, Vhalf, kV, τCa, αCa,
            krel, kclear, kon, koff, kG, kGdeg, G50, nGi,
            kACbasal, kcAMPdeg, kPKA, kPKAdeg) = p

    V,  Ca,  DA,  P,  Gi,  cAMP,  PKA = u   # current state

    # instantaneous Ca current
    m_inf = 1 / (1 + exp(-(V - Vhalf)/kV))
    I_leak = -gL * (V - EL)          # nA (nS·mV)
    I_Ca  = -gCa * m_inf * (V - ECa)        # nA (nS·mV)

    # 1. Voltage
    du[1] = (I_leak + I_Ca + I_stim(t)) / Cm          # dV/dt (mV/ms)

    # 2. Calcium
    du[2] = -Ca/τCa + αCa*I_Ca                               # µM/ms

    # 3. Dopamine
    du[3] = krel*Ca - kclear*DA

    # 4. D2 occupancy
    du[4] = kon*DA*(1 - P) - koff*P

    # 5. Gi/Go
    du[5] = kG*P - kGdeg*Gi

    # 6. cAMP (AC inhibited by Gi)
    AC = kACbasal / (1 + (Gi/G50)^nGi)
    du[6] = AC - kcAMPdeg*cAMP

    # 7. PKA
    du[7] = kPKA*cAMP - kPKAdeg*PKA
end
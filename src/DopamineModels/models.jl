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

"""
    dopaminergic_autoreceptor_grid!(du, u, p, t, grid_params, method)
Modified version of dopaminergic_autoreceptor! that uses a 2D grid for dopamine diffusion.
"""
function dopaminergic_autoreceptor_grid!(du, u, p, t, grid_params, method)
    (Cm, gL, EL, gCa, ECa, Vhalf, kV, τCa, αCa,
            krel, kclear, kon, koff, kG, kGdeg, G50, nGi,
            kACbasal, kcAMPdeg, kPKA, kPKAdeg) = p

    # Extract variables
    V = u[1]
    Ca = u[2]
    DA_grid = u[3:end-4]  # Already a vector
    P = u[end-3]
    Gi = u[end-2]
    cAMP = u[end-1]
    PKA = u[end]

    # instantaneous Ca current
    m_inf = 1 / (1 + exp(-(V - Vhalf)/kV))
    I_leak = -gL * (V - EL)          # nA (nS·mV)
    I_Ca  = -gCa * m_inf * (V - ECa)        # nA (nS·mV)

    # 1. Voltage
    du[1] = (I_leak + I_Ca + I_stim(t)) / Cm          # dV/dt (mV/ms)

    # 2. Calcium
    du[2] = -Ca/τCa + αCa*I_Ca                               # µM/ms

    # 3. Dopamine grid (using the grid discretization)
    update_dopamine_grid!(du[3:end-4], DA_grid, (krel, kclear), Ca, grid_params; method=method)

    # 4. D2 occupancy (using average dopamine concentration)
    DA_avg = mean(DA_grid)
    du[end-3] = kon*DA_avg*(1 - P) - koff*P

    # 5. Gi/Go
    du[end-2] = kG*P - kGdeg*Gi

    # 6. cAMP (AC inhibited by Gi)
    AC = kACbasal / (1 + (Gi/G50)^nGi)
    du[end-1] = AC - kcAMPdeg*cAMP

    # 7. PKA
    du[end] = kPKA*cAMP - kPKAdeg*PKA
end
dRh_dt = (Rh, Rhi) -> Jhm - a1 * Rh + a2 * Rhi
dRhi_dt = (Rh, Rhi) -> a1 * Rh - (a2 + a3) * Rhi
dTr_dt = (Rh, Tr, PDE) -> e * Rh * (Ttot - Tr) - b1 * Tr + s2 * PDE + s1 * Tr * (PDEtot - PDE)
dPDE_dt = (Tr, PDE) -> s1 * Tr * (PDEtot - PDE) - s2 * PDE
dCa_dt = (J, Ca, Cab) -> b * J - cCa * (Ca - C0) - k1 * (eT - Cab) * Ca + k2 * Cab
dCab_dt = (Ca, Cab) -> k1 * (eT - Cab) * Ca - k2 * Cab
dcGMP_dt = (Ca, cGMP) -> Amax / (1 + (Ca / Kc)^4) - cGMP * (V + r * PDE)
Iphoto = (V, cGMP) -> -J_max * (cGMP^3) / (cGMP^3 + 10^3) * (1.0 - exp((V - 8.5) / 17.0))
IKCa = (V, Ca, mKCa) -> gKCa * mKCa^2 * (Ca / (Ca + 0.3)) * (V - EK)

Ih = (V, O1, O2, O3) -> gh * (O1 + O2 + O3) * (V - Eh)
IL = (V) -> gL * (V - EL)
IKv = (V, mKv, hKv) -> gKv * mKv^3 * hKv * (V - EK)
ICa = (V, mCa, Ca_s, Ca_o) -> gCa * mCa^4 * (V - (-12.5 * log10(Ca_s / Ca_o)))
IClCa = (Ca_s, V) -> gCl * (1 / (1 + exp((1 / (0.37 - Ca_s)) / 0.09))) * (V - ECl)

dmKCa_dt = (V, mKCa) -> exp(15 * (-(80 - V) / 40)) * (1 - mKCa) - 20 * exp(-V / 35) * mKCa
dmKv_dt = (V, mKv) -> 0.15 * exp((100 - V) / 42) * (1 - mKv) - 9 * exp(-V / 125) * mKv
dhKv_dt = (V, hKv) -> exp(-V / 120) * (1 - hKv) - 0.1 * exp(-(V + 38) / 7) * hKv
dmCa_dt = (V, mCa) -> 3 * exp((80 - V) / 25) * (1 - mCa) - (1 / (1 + exp((V + 38) / 7))) * mCa

dC1_dt = (V, C1, C2, O1, O2, O3) -> 4 * exp(-(V + 78) / 14) * C2 - (3 * exp(-(V + 78) / 14) + exp(-(V - 8) / 19)) * C1
dC2_dt = (V, C1, C2, O1, O2) -> 3 * exp(-(V + 78) / 14) * C1 - (2 * exp(-(V + 78) / 14) + 2 * exp(-(V - 8) / 19)) * C2 + 2 * exp(-(V - 8) / 19) * O1
dO1_dt = (V, C2, O1, O2) -> 2 * exp(-(V + 78) / 14) * C2 - (exp(-(V + 78) / 14) + 3 * exp(-(V - 8) / 19)) * O1 + 3 * exp(-(V - 8) / 19) * O2
dO2_dt = (V, O1, O2, O3) -> exp(-(V + 78) / 14) * O1 - (4 * exp(-(V - 8) / 19)) * O2 + 4 * exp(-(V - 8) / 19) * O3
dO3_dt = (V, O2, O3) -> exp(-(V + 78) / 14) * O2 - 4 * exp(-(V - 8) / 19) * O3


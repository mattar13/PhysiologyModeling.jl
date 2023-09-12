#Load the auxillary equations
println("Loading auxillary functions")
M∞(v) = (1 + tanh((v - V1) / V2)) / 2
N∞(v) = (1 + tanh((v - V3) / V4)) / 2
Λ(V) = cosh((V - V3) / (2 * V4));
Φe(v) = 1 / (1 + exp(-VSe * (v - V0e)))
Φi(v) = 1 / (1 + exp(-VSi * (v - V0i)))
ħe(e) = (e^2) / (e^2 + k_ACh)
ħi(i) = (i^2) / (i^2 + k_GABA)

α_M(v) = -(v - V8) / (V7 * (exp(-(v - V8) / V9) - 1))
β_M(v) = V10 * (exp(-(v - V11) / V12))
α_H(v) = V13 * (exp(-(v - V14) / V15))
β_H(v) = 1 / (V16 * (exp(-(v - V17) / V18) + 1))

ILeak(v) = -g_leak * (v - E_leak)
ICa(v) = -g_Ca * M∞(v) * (v - E_Ca)
IK(v) = -g_K * n * (v - E_K)
ITREK(v) = -g_TREK * b * (v - E_K)
IACh(v) = -g_ACh * ħe(e) * (v - E_ACh)
IGABA(v) = -g_GABA * ħi(i) * (v - E_Cl)
INa(v) = -g_Na * m^3 * h * (v - E_Na)


# These equations are for the inplace versions
M∞(v, V1, V2) = (1 + tanh((v - V1) / V2)) / 2
N∞(v, V3, V4) = (1 + tanh((v - V3) / V4)) / 2
Λ(V, V3, V4) = cosh((V - V3) / (2 * V4));
Φe(v, VSe, V0e) = 1 / (1 + exp(-VSe * (v - V0e)))
Φi(v, VSi, V0i) = 1 / (1 + exp(-VSi * (v - V0i)))
ħe(e, k_ACh) = (e^2) / (e^2 + k_ACh)
ħi(i, k_GABA) = (i^2) / (i^2 + k_GABA)

α_M(v, V7, V8, V9) = -(v - V8) / (V7 * (exp(-(v - V8) / V9) - 1))
β_M(v, V10, V11, V12) = V10 * (exp(-(v - V11) / V12))
α_H(v, V13, V14, V15) = V13 * (exp(-(v - V14) / V15))
β_H(v, V16, V17, V18) = 1 / (V16 * (exp(-(v - V17) / V18) + 1))

R∞(V, V1, V2) = (1 + tanh((V - V1) / V2)) / 2
fI(g, R, v, E) = -g * R * (v - E)

#ODE==============================================================================================================#

SAC_eqs = [
     Dt(I_ext) ~ -I_ext + I_app, #This parameter is controlled by an outside      
     Dt(v) ~ (ILeak(v) + ICa(v) + IK(v) + ITREK(v) + IACh(v) + IGABA(v) + INa(v) + I_ext + g_W*W) / C_m,
     Dt(n) ~ (Λ(v) * ((N∞(v) - n))) / τn,
     Dt(m) ~ α_M(v) * (1 - m) - β_M(v) * m,
     Dt(h) ~ α_H(v) * (1 - h) - β_H(v) * h,
     Dt(c) ~ (C_0 + δ * (ICa(v)) - λ * c) / τc,
     Dt(a) ~ (α * c^4 * (1 - a) - a) / τa,
     Dt(b) ~ (β * a^4 * (1 - b) - b) / τb,
     Dt(e) ~ (ρe * Φe(v) - e) / τACh, #ρe-e, #
     Dt(i) ~ (ρi * Φi(v) - i) / τGABA, #ρi-i, #
     Dt(W) ~ -W / τw
]

#SDE==============================================================================================================#

SAC_noise_eqs = zeros(length(SAC_eqs))
SAC_noise_eqs[end] = 1.0
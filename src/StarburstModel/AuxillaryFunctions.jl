#Load the auxillary equations
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

ILeak(v, g_leak, E_leak) = -g_leak * (v - E_leak)
ICa(v, g_Ca, V1, V2, E_Ca) = -g_Ca * M∞(v, V1, V2) * (v - E_Ca)
IK(v, n, g_K, E_K) = -g_K * n * (v - E_K)
ITREK(v, b, g_TREK, E_K) = -g_TREK * b * (v - E_K)
IACh(v, e, g_ACh, k_ACh, E_ACh) = -g_ACh * ħe(e, k_ACh) * (v - E_ACh)
IGABA(v, i, g_GABA, k_GABA, E_Cl) = -g_GABA * ħi(i, k_GABA) * (v - E_Cl)
INa(v, m, h, g_Na, E_Na) = -g_Na * m^3 * h * (v - E_Na)

#This is our non-linear distance function
δX(x, a, b, c) = a * exp(-((x - b)^2) / (2 * c^2))
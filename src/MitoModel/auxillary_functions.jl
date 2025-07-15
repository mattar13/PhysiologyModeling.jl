# Auxiliary function for glucose transport
kGLUT(GLU, kMAX; GLUT_MAX = 10.0) = kMAX * (1.0 - GLU / GLUT_MAX)

# Hodgkin-Huxley auxiliary functions
α_m(V) = 0.1 * (V + 40) / (1 - exp(-(V + 40) / 10))
β_m(V) = 4 * exp(-(V + 65) / 18)
α_h(V) = 0.07 * exp(-(V + 65) / 20)
β_h(V) = 1 / (1 + exp(-(V + 35) / 10))
α_n(V) = 0.01 * (V + 55) / (1 - exp(-(V + 55) / 10))
β_n(V) = 0.125 * exp(-(V + 65) / 80)

m_inf(V) = α_m(V) / (α_m(V) + β_m(V))
h_inf(V) = α_h(V) / (α_h(V) + β_h(V))
n_inf(V) = α_n(V) / (α_n(V) + β_n(V))
τ_m(V) = 1 / (α_m(V) + β_m(V))
τ_h(V) = 1 / (α_h(V) + β_h(V))
τ_n(V) = 1 / (α_n(V) + β_n(V))

# Ionic currents
I_Na(V, m, h, g_Na, E_Na) = g_Na * m^3 * h * (V - E_Na)
I_K(V, n, g_K, E_K) = g_K * n^4 * (V - E_K)
I_leak(V, g_leak, E_leak) = g_leak * (V - E_leak)

# ATP-dependent functions
g_Na_ATP(ATP, g_Na_max, K_ATP_Na) = g_Na_max * ATP / (K_ATP_Na + ATP)
g_K_ATP(ATP, g_K_max, K_ATP_K) = g_K_max * ATP / (K_ATP_K + ATP)
I_NaK_ATPase(ATP, I_NaK_max, K_ATP_NaK) = I_NaK_max * ATP / (K_ATP_NaK + ATP)

# Time-dependent stimulus function
I_APP(t, stim_start, stim_end, I_amplitude) = (t >= stim_start && t <= stim_end) ? I_amplitude : 0.0

# Register functions with Catalyst
@register_symbolic g_Na_ATP(ATP, g_Na_max, K_ATP_Na)
@register_symbolic g_K_ATP(ATP, g_K_max, K_ATP_K)
@register_symbolic I_NaK_ATPase(ATP, I_NaK_max, K_ATP_NaK)
@register_symbolic I_APP(t, stim_start, stim_end, I_amplitude)

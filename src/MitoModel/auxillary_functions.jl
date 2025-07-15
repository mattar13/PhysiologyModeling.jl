# Auxiliary function for glucose transport
kGLUT(GLU, kMAX, GLUT_MAX) = kMAX * (1.0 - GLU / GLUT_MAX)
kOXYGEN(O2, kMAX, O2_MAX) = kMAX * (1.0 - O2 / O2_MAX)

# Hodgkin-Huxley auxiliary functions
alpha_m(V) = 0.1 * (V + 40) / (1 - exp(-(V + 40) / 10))
beta_m(V) = 4 * exp(-(V + 65) / 18)
alpha_h(V) = 0.07 * exp(-(V + 65) / 20)
beta_h(V) = 1 / (1 + exp(-(V + 35) / 10))
alpha_n(V) = 0.01 * (V + 55) / (1 - exp(-(V + 55) / 10))
beta_n(V) = 0.125 * exp(-(V + 65) / 80)

# Ionic currents
I_Na(V, m, h, g_Na, E_Na) = g_Na * m^3 * h * (V - E_Na)
I_K(V, n, g_K, E_K) = g_K * n^4 * (V - E_K)
I_leak(V, g_leak, E_leak) = g_leak * (V - E_leak)

# Time-dependent stimulus function
I_APP(t, stim_start, stim_end, I_amplitude) = (t >= stim_start && t <= stim_end) ? I_amplitude : 0.0

#Function for ATP breakdown due to the NaK pump
h_PUMP(ATP, K_ATP_NaK) = ATP / (K_ATP_NaK + ATP)
k_PUMP(V) = 1 + 0.5 * (1 / (1 + exp(-(V + 40) / 10)))  # Enhanced at depolarized V
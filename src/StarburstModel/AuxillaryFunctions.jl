#Load the auxillary equations
@inline M∞(v::T, V1::Float64, V2::Float64) where T = (1 + tanh((v - V1) / V2)) / 2
@inline N∞(v::T, V3::Float64, V4::Float64) where T = (1 + tanh((v - V3) / V4)) / 2
@inline Λ(V::T, V3::Float64, V4::Float64) where T = cosh((V - V3) / (2 * V4));
@inline Φe(v::T, VSe::Float64, V0e::Float64) where T = 1 / (1 + exp(-VSe * (v - V0e)))
@inline Φi(v::T, VSi::Float64, V0i::Float64) where T = 1 / (1 + exp(-VSi * (v - V0i)))
@inline ħe(e::T, k_ACh::Float64) where T = (e^2) / (e^2 + k_ACh)
@inline ħi(i::T, k_GABA::Float64) where T = (i^2) / (i^2 + k_GABA)
@inline ħg(g::T, k_GLUT::Float64) where T = (g^2) / (g^2 + k_GLUT)

@inline α_M(v::T, V7::Float64, V8::Float64, V9::Float64) where T = -(v - V8) / (V7 * (exp(-(v - V8) / V9) - 1))
@inline β_M(v::T, V10::Float64, V11::Float64, V12::Float64) where T = V10 * (exp(-(v - V11) / V12))
@inline α_H(v::T, V13::Float64, V14::Float64, V15::Float64) where T = V13 * (exp(-(v - V14) / V15))
@inline β_H(v::T, V16::Float64, V17::Float64, V18::Float64) where T = 1 / (V16 * (exp(-(v - V17) / V18) + 1))

@inline ILeak(v::T, g_leak::Float64, E_leak::Float64) where T = -g_leak * (v - E_leak)
@inline ICa(v::T, g_Ca::Float64, V1::Float64, V2::Float64, E_Ca::Float64) where T = -g_Ca * M∞(v, V1, V2) * (v - E_Ca)
@inline IK(v::T, n::T, g_K::Float64, E_K::Float64) where T = -g_K * n * (v - E_K)
@inline ITREK(v::T, b::T, g_TREK::Float64, E_K::Float64) where T = -g_TREK * b * (v - E_K)
@inline IACh(v::T, e::T, g_ACh::Float64, k_ACh::Float64, E_ACh::Float64) where T = -g_ACh * ħe(e, k_ACh) * (v - E_ACh)
@inline IGABA(v::T, i::T, g_GABA::Float64, k_GABA::Float64, E_Cl::Float64) where T = -g_GABA * ħi(i, k_GABA) * (v - E_Cl)
@inline INa(v::T, m::T, h::T, g_Na::Float64, E_Na::Float64) where T = -g_Na * m^3 * h * (v - E_Na)

#==================================================[Equations for Glutamate dynamics]==================================================#
@inline ICa_mGluR2(v::T, q::T, g_Ca::Float64, V1::Float64, V2::Float64, E_Ca::Float64) where T = -g_Ca * M∞(v, V1, V2) * (1.0-q) * (v - E_Ca)
@inline IGLUT(v::T, g::T, g_GLUT::Float64, k_GLUT::Float64, E_GLUT::Float64) where T = -g_GLUT * ħe(g, k_GLUT) * (v - E_GLUT)

gauss_pulse(t; t0::Float64  = 25.0, spread::Float64 = 500.0, peak_amp::Float64 = 1.0) = peak_amp * exp(-(t0-t)^2/((2*spread)^2))
#Load the auxillary equations
@inline M∞(v::T, V1, V2) where T = (1 + tanh((v - V1) / V2)) / 2
@inline N∞(v::T, V3, V4) where T = (1 + tanh((v - V3) / V4)) / 2
@inline Λ(V::T, V3, V4) where T = cosh((V - V3) / (2 * V4));
@inline Φe(v::T, VSe, V0e) where T = 1 / (1 + exp(-VSe * (v - V0e)))
@inline Φi(v::T, VSi, V0i) where T = 1 / (1 + exp(-VSi * (v - V0i)))
@inline ħe(e::T, k_ACh) where T = (e^2) / (e^2 + k_ACh)
@inline ħi(i::T, k_GABA) where T = (i^2) / (i^2 + k_GABA)
@inline ħg(g::T, k_GLUT) where T = (g^2) / (g^2 + k_GLUT)

@inline α_M(v::T, V7, V8, V9) where T = -(v - V8) / (V7 * (exp(-(v - V8) / V9) - 1))
@inline β_M(v::T, V10, V11, V12) where T = V10 * (exp(-(v - V11) / V12))
@inline α_H(v::T, V13, V14, V15) where T = V13 * (exp(-(v - V14) / V15))
@inline β_H(v::T, V16, V17, V18) where T = 1 / (V16 * (exp(-(v - V17) / V18) + 1))

@inline ILeak(v::T, g_leak, E_leak) where T = -g_leak * (v - E_leak)
@inline ICa(v::T, g_Ca, V1, V2, E_Ca) where T = -g_Ca * M∞(v, V1, V2) * (v - E_Ca)
#This function is for the potassium current. g_K can be a vector. 
@inline IK(v::T, n::T, g_K::T, E_K) where T = -g_K .* n .* (v .- E_K)
@inline IK(v::T, n::T, g_K::eltype(T), E_K) where T = -g_K .* n .* (v .- E_K)

@inline ITREK(v::T, b::T, g_TREK, E_K) where T = -g_TREK * b * (v - E_K)
@inline IACh(v::T, e::T, g_ACh, k_ACh, E_ACh) where T = -g_ACh * ħe(e, k_ACh) * (v - E_ACh)
@inline IGABA(v::T, i::T, g_GABA, k_GABA, E_Cl) where T = -g_GABA * ħi(i, k_GABA) * (v - E_Cl)
@inline INa(v::T, m::T, h::T, g_Na, E_Na) where T = -g_Na * m^3 * h * (v - E_Na)

#==================================================[Equations for Glutamate dynamics]==================================================#
@inline ICa_mGluR2(v::T, q::T, g_Ca, V1, V2, E_Ca) where T = -g_Ca * M∞(v, V1, V2) * (1.0-q) * (v - E_Ca)
@inline IGLUT(v::T, g::T, g_GLUT, k_GLUT, E_GLUT) where T = -g_GLUT * ħe(g, k_GLUT) * (v - E_GLUT)

gauss_pulse(t; t0  = 25.0, spread = 500.0, peak_amp = 1.0) = peak_amp * exp(-(t0-t)^2/((2*spread)^2))

#%% Callback functions 
function apply_glutamate_affect!(integrator; xs = nothing, ys = nothing,  
    pulse_list = [93, 45, 21, 9, 3, 1, 2, 6, 12, 30, 63], 
    n_stops = 10,x_stops = nothing,
    dt_pulse = 250.0, pulse_start = 200.0, 
    spread = 250.0, amp = 5.0
) 
    if !isnothing(x_stops)
        for i in 1:n_stops-1
            for idx in findall(x_stops[i] .<= xs .<= x_stops[i+1])
                integrator.u[idx, 11] = gauss_pulse(integrator.t; 
                    t0 = pulse_start + (dt_pulse * i), 
                    spread = spread, peak_amp = amp)
            end
        end
    else
        for (i, pulse) in enumerate(pulse_list)
            integrator.u[pulse, 11] = gauss_pulse(integrator.t; t0 = pulse_start + (dt_pulse * i), spread = spread, peak_amp = amp)
        end
    end
end
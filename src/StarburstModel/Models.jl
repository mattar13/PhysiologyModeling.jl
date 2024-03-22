#=================================================ODE EQUATIONS (Basics)=================================================#
function SAC_ODE(du, u, p, t)
     I_ext = view(u, 1)
     v = view(u, 2)
     n = view(u, 3)
     m = view(u, 4)
     h = view(u, 5)
     c = view(u, 6)
     a = view(u, 7)
     b = view(u, 8)
     e = view(u, 9)
     i = view(u, 10)
     g = view(u, 11)
     q = view(u, 12)
     W = view(u, 13)

     dI_ext = view(du, 1)
     dv = view(du, 2)
     dn = view(du, 3)
     dm = view(du, 4)
     dh = view(du, 5)
     dc = view(du, 6)
     da = view(du, 7)
     db = view(du, 8)
     de = view(du, 9)
     di = view(du, 10)
     dg = view(du, 11)
     dq = view(du, 12)
     dW = view(du, 13)

     (I_app, VC,
          C_m, g_W, τw, 
          g_leak, E_leak, 
          g_K, V3, V4, E_K, τn, 
          g_Ca, V1, V2,E_Ca, τc,
          g_Na, E_Na, 
          g_TREK,
          C_0, λ , δ,  
          α, τa, 
          β, τb, 
          a_n, b_n,
          VSe, V0e, ρe,  g_ACh, k_ACh, E_ACh,  τACh,
          VSi, V0i, ρi,  g_GABA, k_GABA, E_Cl, τGABA, 
          g_GLUT, k_GLUT, E_GLUT, 
          γg, g_n, τq,
          V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18, 
          stim_start, stim_stop
     ) = p

     @. dI_ext = I_app-I_ext
     @. dv = (ILeak(v, g_leak, E_leak) + 
          + ICa_mGluR2(v, q, g_Ca, V1, V2, E_Ca) + IK(v, n, g_K, E_K) + INa(v, m, h, g_Na, E_Na)
          + ITREK(v, b, g_TREK, E_K) 
          + IACh(v, e, g_ACh, k_ACh, E_ACh) 
          + IGABA(v, i, g_GABA, k_GABA, E_Cl) 
          + IGLUT(v, g, g_GLUT, k_GLUT, E_GLUT) #These are ionic glutamate channels
          + I_app + W) / C_m #Unless we are doing IC, this has to stay this way
     @. dn = (Λ(v, V3, V4) * ((N∞(v, V3, V4) - n))) / τn
     @. dm = α_M(v, V7, V8, V9) * (1 - m) - β_M(v, V10, V11, V12) * m
     @. dh = α_H(v, V13, V14, V15) * (1 - h) - β_H(v, V16, V17, V18) * h
     @. dc = (C_0 + δ * (ICa_mGluR2(v, q, g_Ca, V1, V2, E_Ca)) - λ * c) / τc
     @. da = (α * c^a_n * (1 - a) - a) / τa #These were the old options
     @. db = (β * a^b_n * (1 - b) - b) / τb #These were the old options
     @. de = (ρe * Φe(v, VSe, V0e) - e) / τACh
     @. di = (ρi * Φi(v, VSi, V0i) - i) / τGABA
     @. dg = 0.0
     @. dq = (γg*g^g_n * (1-q) - q) / τq
     @. dW = -W / τw
     nothing
end

function SAC_ODE_IC(du, u, p, t; stim_start = 500.0, stim_stop = 2000.0)
     I_ext = view(u, 1)
     v = view(u, 2)
     n = view(u, 3)
     m = view(u, 4)
     h = view(u, 5)
     c = view(u, 6)
     a = view(u, 7)
     b = view(u, 8)
     e = view(u, 9)
     i = view(u, 10)
     g = view(u, 11)
     q = view(u, 12)
     W = view(u, 13)

     dI_ext = view(du, 1)
     dv = view(du, 2)
     dn = view(du, 3)
     dm = view(du, 4)
     dh = view(du, 5)
     dc = view(du, 6)
     da = view(du, 7)
     db = view(du, 8)
     de = view(du, 9)
     di = view(du, 10)
     dg = view(du, 11)
     dq = view(du, 12)
     dW = view(du, 13)

     (I_app, VC,
          C_m, g_W, τw, 
          g_leak, E_leak, 
          g_K, V3, V4, E_K, τn, 
          g_Ca, V1, V2,E_Ca, τc,
          g_Na, E_Na, 
          g_TREK,
          C_0, λ , δ,  
          α, τa, 
          β, τb, 
          a_n, b_n,
          VSe, V0e, ρe,  g_ACh, k_ACh, E_ACh,  τACh,
          VSi, V0i, ρi,  g_GABA, k_GABA, E_Cl, τGABA, 
          g_GLUT, k_GLUT, E_GLUT, 
          γg, g_n, τq,
          V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18, 
          stim_start, stim_stop
     ) = p
     
     if stim_start < t < stim_stop
          stim_amp = I_app
     else
          stim_amp = 0.0
     end

     @. dI_ext = (stim_amp-I_ext)/0.001
     @. dv = (ILeak(v, g_leak, E_leak) + 
          + ICa_mGluR2(v, q, g_Ca, V1, V2, E_Ca) + IK(v, n, g_K, E_K) + INa(v, m, h, g_Na, E_Na)
          + ITREK(v, b, g_TREK, E_K) 
          + IACh(v, e, g_ACh, k_ACh, E_ACh) 
          + IGABA(v, i, g_GABA, k_GABA, E_Cl) 
          + IGLUT(v, g, g_GLUT, k_GLUT, E_GLUT) #These are ionic glutamate channels
          + stim_amp + W) / C_m #Unless we are doing IC, this has to stay this way
     @. dn = (Λ(v, V3, V4) * ((N∞(v, V3, V4) - n))) / τn
     @. dm = α_M(v, V7, V8, V9) * (1 - m) - β_M(v, V10, V11, V12) * m
     @. dh = α_H(v, V13, V14, V15) * (1 - h) - β_H(v, V16, V17, V18) * h
     @. dc = (C_0 + δ * (ICa_mGluR2(v, q, g_Ca, V1, V2, E_Ca)) - λ * c) / τc
     @. da = (α * c^a_n * (1 - a) - a) / τa #These were the old options
     @. db = (β * a^b_n * (1 - b) - b) / τb #These were the old options
     @. de = (ρe * Φe(v, VSe, V0e) - e) / τACh
     @. di = (ρi * Φi(v, VSi, V0i) - i) / τGABA
     @. dg = 0.0
     @. dq = (γg*g^g_n * (1-q) - q) / τq
     @. dW = -W / τw
     nothing
end

function SAC_ODE_VC(du, u, p, t; stim_start = 500.0, stim_stop = 2000.0, hold = -65.0, k = 1000.0)
     I_ext = view(u, 1)
     v = view(u, 2)
     n = view(u, 3)
     m = view(u, 4)
     h = view(u, 5)
     c = view(u, 6)
     a = view(u, 7)
     b = view(u, 8)
     e = view(u, 9)
     i = view(u, 10)
     g = view(u, 11)
     q = view(u, 12)
     W = view(u, 13)

     dI_ext = view(du, 1)
     dv = view(du, 2)
     dn = view(du, 3)
     dm = view(du, 4)
     dh = view(du, 5)
     dc = view(du, 6)
     da = view(du, 7)
     db = view(du, 8)
     de = view(du, 9)
     di = view(du, 10)
     dg = view(du, 11)
     dq = view(du, 12)
     dW = view(du, 13)

     (I_app, VC,
          C_m, g_W, τw, 
          g_leak, E_leak, 
          g_K, V3, V4, E_K, τn, 
          g_Ca, V1, V2,E_Ca, τc,
          g_Na, E_Na, 
          g_TREK,
          C_0, λ , δ,  
          α, τa, 
          β, τb, 
          a_n, b_n,
          VSe, V0e, ρe,  g_ACh, k_ACh, E_ACh,  τACh,
          VSi, V0i, ρi,  g_GABA, k_GABA, E_Cl, τGABA, 
          g_GLUT, k_GLUT, E_GLUT, 
          γg, g_n, τq,
          V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18, 
          stim_start, stim_stop
     ) = p
     
     if stim_start < t < stim_stop
          @. dv = (VC-v)*k
     else
          @. dv = (hold-v)*k
     end
     @. dI_ext = (ILeak(v, g_leak, E_leak) + 
          + ICa_mGluR2(v, q, g_Ca, V1, V2, E_Ca) + IK(v, n, g_K, E_K) + INa(v, m, h, g_Na, E_Na)
          + ITREK(v, b, g_TREK, E_K) 
          + IACh(v, e, g_ACh, k_ACh, E_ACh) 
          + IGABA(v, i, g_GABA, k_GABA, E_Cl) 
          + IGLUT(v, g, g_GLUT, k_GLUT, E_GLUT) #These are ionic glutamate channels
     ) - I_ext
     @. dn = (Λ(v, V3, V4) * ((N∞(v, V3, V4) - n))) / τn
     @. dm = α_M(v, V7, V8, V9) * (1 - m) - β_M(v, V10, V11, V12) * m
     @. dh = α_H(v, V13, V14, V15) * (1 - h) - β_H(v, V16, V17, V18) * h
     @. dc = (C_0 + δ * (ICa_mGluR2(v, q, g_Ca, V1, V2, E_Ca)) - λ * c) / τc
     @. da = (α * c^a_n * (1 - a) - a) / τa #These were the old options
     @. db = (β * a^b_n * (1 - b) - b) / τb #These were the old options
     @. de = (ρe * Φe(v, VSe, V0e) - e) / τACh
     @. di = (ρi * Φi(v, VSi, V0i) - i) / τGABA
     @. dg = 0.0
     @. dq = (γg*g^g_n * (1-q) - q) / τq
     @. dW = -W / τw
     nothing
end

#=================================================PDE EQUATIONS=================================================#
function ∇α(du, u, cell_map, t) #Could it really be this easy? 
     du .+= (cell_map.strength_out .* u) + (cell_map.strength * u)
end

function DIFFUSION_MODEL(du, u, p, t; active_cell = 221, growth_rate = 0.5)
     du .= -u/540 #du decays over time
     if 500.0 < t < 2500.0
          #We want to add some diffusive material during a time range
          du[active_cell] = growth_rate
     end
     ∇α(du, u, p, t)#Diffusion occurs after
     #We should go through and decay the edges 
end

function SAC_PDE(du, u, p, t, MAP) 
     I_ext = view(u, :, 1)
     v = view(u, :, 2)
     n = view(u, :, 3)
     m = view(u, :, 4)
     h = view(u, :, 5)
     c = view(u, :, 6)
     a = view(u, :, 7)
     b = view(u, :, 8)
     e = view(u, :, 9)
     i = view(u, :, 10)
     g = view(u, :, 11)
     q = view(u, :, 12)
     W = view(u, :, 13)

     dI_ext = view(du, :, 1)
     dv = view(du, :, 2)
     dn = view(du, :, 3)
     dm = view(du, :, 4)
     dh = view(du, :, 5)
     dc = view(du, :, 6)
     da = view(du, :, 7)
     db = view(du, :, 8)
     de = view(du, :, 9)
     di = view(du, :, 10)
     dg = view(du, :, 11)
     dq = view(du, :, 12)
     dW = view(du, :, 13)

     (I_app, VC,
          C_m, g_W, τw, 
          g_leak, E_leak, 
          g_K, V3, V4, E_K, τn, 
          g_Ca, V1, V2,E_Ca, τc,
          g_Na, E_Na, 
          g_TREK,
          C_0, λ , δ,  
          α, τa, 
          β, τb, 
          a_n, b_n,
          VSe, V0e, ρe,  g_ACh, k_ACh, E_ACh,  τACh,
          VSi, V0i, ρi,  g_GABA, k_GABA, E_Cl, τGABA, 
          g_GLUT, k_GLUT, E_GLUT, 
          γg, g_n, τq,
          V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18, 
          stim_start, stim_stop
     ) = p

     @. dI_ext = I_app-I_ext
     @. dv = (ILeak(v, g_leak, E_leak) + 
          + ICa_mGluR2(v, q, g_Ca, V1, V2, E_Ca) + IK(v, n, g_K, E_K) + INa(v, m, h, g_Na, E_Na)
          + ITREK(v, b, g_TREK, E_K) 
          + IACh(v, e, g_ACh, k_ACh, E_ACh) 
          + IGABA(v, i, g_GABA, k_GABA, E_Cl) 
          + IGLUT(v, g, g_GLUT, k_GLUT, E_GLUT) #These are ionic glutamate channels
          + I_app + W) / C_m #Unless we are doing IC, this has to stay this way
     @. dn = (Λ(v, V3, V4) * ((N∞(v, V3, V4) - n))) / τn
     @. dm = α_M(v, V7, V8, V9) * (1 - m) - β_M(v, V10, V11, V12) * m
     @. dh = α_H(v, V13, V14, V15) * (1 - h) - β_H(v, V16, V17, V18) * h
     @. dc = (C_0 + δ * (ICa_mGluR2(v, q, g_Ca, V1, V2, E_Ca)) - λ * c) / τc
     @. da = (α * c^a_n * (1 - a) - a) / τa #These were the old options
     @. db = (β * a^b_n * (1 - b) - b) / τb #These were the old options
     @. de = (ρe * Φe(v, VSe, V0e) - e) / τACh
     ∇α(de, e, MAP, t)
     @. di = (ρi * Φi(v, VSi, V0i) - i) / τGABA
     @. dg = 0.0
     @. dq = (γg*g^g_n * (1-q) - q) / τq
     @. dW = -W / τw
     nothing
end

function SAC_GAP(du, u, p, t, MAP; gGAP = 0.1)
     I_ext = view(u, :, 1)
     v = view(u, :, 2)
     n = view(u, :, 3)
     m = view(u, :, 4)
     h = view(u, :, 5)
     c = view(u, :, 6)
     a = view(u, :, 7)
     b = view(u, :, 8)
     e = view(u, :, 9)
     i = view(u, :, 10)
     g = view(u, :, 11)
     q = view(u, :, 12)
     W = view(u, :, 13)

     dI_ext = view(du, :, 1)
     dv = view(du, :, 2)
     dn = view(du, :, 3)
     dm = view(du, :, 4)
     dh = view(du, :, 5)
     dc = view(du, :, 6)
     da = view(du, :, 7)
     db = view(du, :, 8)
     de = view(du, :, 9)
     di = view(du, :, 10)
     dg = view(du, :, 11)
     dq = view(du, :, 12)
     dW = view(du, :, 13)

     (I_app, VC,
          C_m, g_W, τw, 
          g_leak, E_leak, 
          g_K, V3, V4, E_K, τn, 
          g_Ca, V1, V2,E_Ca, τc,
          g_Na, E_Na, 
          g_TREK,
          C_0, λ , δ,  
          α, τa, 
          β, τb, 
          a_n, b_n,
          VSe, V0e, ρe,  g_ACh, k_ACh, E_ACh,  τACh,
          VSi, V0i, ρi,  g_GABA, k_GABA, E_Cl, τGABA, 
          g_GLUT, k_GLUT, E_GLUT, 
          γg, g_n, τq,
          V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18, 
          stim_start, stim_stop
     ) = p

     @. dI_ext = I_app-I_ext #I am going to selfishly use this. Can figure this out later
     @. dv = (ILeak(v, g_leak, E_leak) + 
          + ICa_mGluR2(v, q, g_Ca, V1, V2, E_Ca) + IK(v, n, g_K, E_K) + INa(v, m, h, g_Na, E_Na)
          + ITREK(v, b, g_TREK, E_K) 
          + IACh(v, e, g_ACh, k_ACh, E_ACh) 
          + IGABA(v, i, g_GABA, k_GABA, E_Cl) 
          + IGLUT(v, g, g_GLUT, k_GLUT, E_GLUT) #These are ionic glutamate channels
          + I_app + W) / C_m #Unless we are doing IC, this has to stay this way
     ∇α(dv, v, MAP, t)

     @. dn = (Λ(v, V3, V4) * ((N∞(v, V3, V4) - n))) / τn
     @. dm = α_M(v, V7, V8, V9) * (1 - m) - β_M(v, V10, V11, V12) * m
     @. dh = α_H(v, V13, V14, V15) * (1 - h) - β_H(v, V16, V17, V18) * h
     @. dc = (C_0 + δ * (ICa_mGluR2(v, q, g_Ca, V1, V2, E_Ca)) - λ * c) / τc
     @. da = (α * c^a_n * (1 - a) - a) / τa #These were the old options
     @. db = (β * a^b_n * (1 - b) - b) / τb #These were the old options
     @. de = (ρe * Φe(v, VSe, V0e) - e) / τACh
     @. di = (ρi * Φi(v, VSi, V0i) - i) / τGABA
     @. dg = 0.0
     @. dq = (γg*g^g_n * (1-q) - q) / τq
     @. dW = -W / τw
     nothing
end

#=================================================NOISE EQUATIONS=================================================#

DIFFUSION_NOISE(du, u, p, t) = du[:] .= 0.001

noise2D(du, u, p, t) = du[:, end] .= p[4]
function noise1D(du, u, p, t) 
     du[end] = p[4]
end

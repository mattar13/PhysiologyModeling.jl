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
     d = view(u, 12)
     q = view(u, 13)
     W = view(u, 14)

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
     dd = view(du, 12)
     dq = view(du, 13)
     dW = view(du, 14)

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
          k_SYT, n_SYT,
          ρe,  g_ACh, k_ACh, E_ACh,  τACh,
          ρi,  g_GABA, k_GABA, E_Cl, τGABA, 
          g_GLUT, k_GLUT, E_GLUT, 
          τd,
          γg, g_n, τq,
          V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18, 
          stim_start, stim_stop
     ) = p
     @. dI_ext = (I_app-I_ext)
     @. dv = (ILeak(v, g_leak, E_leak) + 
          + ICa(v, g_Ca, V1, V2, E_Ca)# * (1.0-q)
          + IK(v, n, g_K, E_K) + INa(v, m, h, g_Na, E_Na)
          + ITREK(v, b, g_TREK, E_K) 
          + IACh(v, e, g_ACh, k_ACh, E_ACh) 
          + IGABA(v, i, g_GABA, k_GABA, E_Cl) 
          + IGLUT(v, g, g_GLUT, k_GLUT, E_GLUT) #These are ionic glutamate channels
          + I_ext + W) / C_m #Unless we are doing IC, this has to stay this way
     @. dn = (Λ(v, V3, V4) * ((N∞(v, V3, V4) - n))) / τn
     @. dm = α_M(v, V7, V8, V9) * (1 - m) - β_M(v, V10, V11, V12) * m
     @. dh = α_H(v, V13, V14, V15) * (1 - h) - β_H(v, V16, V17, V18) * h
     @. dc = (C_0 + δ * (ICa(v, g_Ca, V1, V2, E_Ca)) - λ * c) / τc
     
     #@. da = (α * c^a_n * (1 - a) - a) / τa #These were the old options
     #@. db = (β * a^b_n * (1 - b) - b) / τb #These were the old options
     @. da = (-α*c^a_n*a + (1-a)) / τa #Why didn't the old way work? 
     @. db = (β*(1-a)^b_n*(1-b) - b)/τb #TREK activity

     #What if we make these dependent on the calcium concentration 
     @. de = (ρe * ΦCa(c, k_SYT, n_SYT) - e) / τACh #Change these to reflect calcium
     @. di = (ρi * ΦCa(c, k_SYT, n_SYT) - i) / τGABA
     @. dg = 0.0 #This value needs to exponentially decay
     @. dd = -d/τd #This is the reuptake and removal of dopamine from synapses
     @. dq = (γg*d^g_n * (1-q) - q) / τq #This is the DA
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
          k_SYT, n_SYT,
          ρe,  g_ACh, k_ACh, E_ACh,  τACh,
          ρi,  g_GABA, k_GABA, E_Cl, τGABA, 
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

     @. dI_ext = (stim_amp-I_ext)
     @. dv = (ILeak(v, g_leak, E_leak) + 
          + ICa(v, g_Ca, V1, V2, E_Ca) * (1.0-q)
          + IK(v, n, g_K, E_K) + INa(v, m, h, g_Na, E_Na)
          + ITREK(v, b, g_TREK, E_K) 
          + IACh(v, e, g_ACh, k_ACh, E_ACh) 
          + IGABA(v, i, g_GABA, k_GABA, E_Cl) 
          + IGLUT(v, g, g_GLUT, k_GLUT, E_GLUT) #These are ionic glutamate channels
          + I_ext + W) / C_m #Unless we are doing IC, this has to stay this way
     @. dn = (Λ(v, V3, V4) * ((N∞(v, V3, V4) - n))) / τn
     @. dm = α_M(v, V7, V8, V9) * (1 - m) - β_M(v, V10, V11, V12) * m
     @. dh = α_H(v, V13, V14, V15) * (1 - h) - β_H(v, V16, V17, V18) * h
     @. dc = (C_0 + δ * (ICa(v,  g_Ca, V1, V2, E_Ca)*(1.0-q)) - λ * c) / τc
     @. da = (α * c^a_n * (1 - a) - a) / τa #These were the old options
     @. db = (β * a^b_n * (1 - b) - b) / τb #These were the old options
     @. de = (ρe * ΦCa(c, k_SYT, n_SYT) - e) / τACh #Change these to reflect calcium
     @. di = (ρi * ΦCa(c, k_SYT, n_SYT) - i) / τGABA
     @. dg = 0.0
     @. dq = (γg*g^g_n * (1-q) - q) / τq
     @. dW = -W / τw
     nothing
end

function SAC_ODE_VC(du, u, p, t; hold = -65.0, k = 1.0)
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
      
     @. dv = (ILeak(v, g_leak, E_leak) + 
          + ICa(v, g_Ca, V1, V2, E_Ca) * (1.0-q)
          + INa(v, m, h, g_Na, E_Na)
          + IK(v, n, g_K, E_K) 
          + ITREK(v, b, g_TREK, E_K) 
          + IACh(v, e, g_ACh, k_ACh, E_ACh) 
          + IGABA(v, i, g_GABA, k_GABA, E_Cl) 
          + IGLUT(v, g, g_GLUT, k_GLUT, E_GLUT) #These are ionic glutamate channels
          + I_app - I_ext + W
     ) / C_m #Unless we are doing IC, this has to stay this way

     if stim_start < t < stim_stop
          @. dI_ext = k*(v+dv - VC) - I_ext
     else
          @. dI_ext = k*(v+dv - hold) - I_ext 
     end   

     @. dn = (Λ(v, V3, V4) * ((N∞(v, V3, V4) - n))) / τn
     @. dm = α_M(v, V7, V8, V9) * (1 - m) - β_M(v, V10, V11, V12) * m
     @. dh = α_H(v, V13, V14, V15) * (1 - h) - β_H(v, V16, V17, V18) * h
     @. dc = (C_0 + δ * (ICa(v,  g_Ca, V1, V2, E_Ca)*(1.0-q)) - λ * c) / τc
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
          u[active_cell] = 1.0
     end
     ∇α(du, u, p, t)#Diffusion occurs after
     #We should go through and decay the edges 
end

function SAC_PDE(du, u, p, t, G_MAP, E_MAP, I_MAP) 
     v_GAP = view(u, :, 1)
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
     d = view(u, :, 13)
     W = view(u, :, 14)

     dv_GAP = view(du, :, 1)
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
     dd = view(du, :, 13)
     dW = view(du, :, 14)

     (I_app, VC, g_GAP,
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
          k_SYT, n_SYT,
          ρe,  g_ACh, k_ACh, E_ACh,  τACh,
          ρi,  g_GABA, k_GABA, E_Cl, τGABA, 
          g_GLUT, k_GLUT, E_GLUT, 
          τd,
          γg, g_n, τq,
          V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18, 
          stim_start, stim_stop
     ) = p

     @. dv_GAP = -v_GAP #Lets try this. 
     ∇α(dv_GAP, v, G_MAP, t)

     @. dv = (ILeak(v, g_leak, E_leak) + 
          + ICa(v, g_Ca, V1, V2, E_Ca) * (1.0-q)
          + IK(v, n, g_K, E_K) + INa(v, m, h, g_Na, E_Na)
          + ITREK(v, b, g_TREK, E_K) 
          + IACh(v, e, g_ACh, k_ACh, E_ACh) 
          + IGABA(v, i, g_GABA, k_GABA, E_Cl) 
          + IGLUT(v, g, g_GLUT, k_GLUT, E_GLUT) #These are ionic glutamate channels
          + I_app 
          + -g_GAP*(1.0-a^2)*v_GAP 
          + W) / C_m #Unless we are doing IC, this has to stay this way
     @. dn = (Λ(v, V3, V4) * ((N∞(v, V3, V4) - n))) / τn
     @. dm = α_M(v, V7, V8, V9) * (1 - m) - β_M(v, V10, V11, V12) * m
     @. dh = α_H(v, V13, V14, V15) * (1 - h) - β_H(v, V16, V17, V18) * h
     @. dc = (C_0 + δ * (ICa(v,  g_Ca, V1, V2, E_Ca)*(1.0-q)) - λ * c) / τc
     
     #@. da = (α * c^a_n * (1 - a) - a) / τa #These were the old options
     #@. db = (β * a^b_n * (1 - b) - b) / τb #These were the old options
     @. da = (-α*c^a_n*a + (1-a)) / τa #Why didn't the old way work? 
     @. db = (β*(1-a)^b_n*(1-b) - b)/τb #TREK activity

     @. de = (ρe * ΦCa(c, k_SYT, n_SYT) - e) / τACh #Change these to reflect calcium
     ∇α(de, e, E_MAP, t)
     @. di = (ρi * ΦCa(c, k_SYT, n_SYT) - i) / τGABA
     ∇α(di, i, I_MAP, t)

     @. dg = -g/500.0
     @. dd = -d / τd
     @. dq = (γg*g^g_n * (1-q) - q) / τq
     @. dW = -W / τw
     nothing
end

#=================================================NOISE EQUATIONS=================================================#

DIFFUSION_NOISE(du, u, p, t) = du[:] .= 0.001

noise2D(du, u, p, t) = du[:, end] .= p[5]
function noise1D(du, u, p, t) 
     du[end] = p[5]
end


"""
We have to remove the two variables that don't change
"""
function DynamicSAC(du, u, p, t)
     #Try making functions for differential equations
     v = view(u, 1)
     n = view(u, 2)
     m = view(u, 3)
     h = view(u, 4)
     c = view(u, 5)
     a = view(u, 6)
     b = view(u, 7)

     dv = view(du, 1)
     dn = view(du, 2)
     dm = view(du, 3)
     dh = view(du, 4)
     dc = view(du, 5)
     da = view(du, 6)
     db = view(du, 7)

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

     @. dv = (ILeak(v, g_leak, E_leak) + 
          + ICa(v, g_Ca, V1, V2, E_Ca)
          + IK(v, n, g_K, E_K) + INa(v, m, h, g_Na, E_Na)
          + ITREK(v, b, g_TREK, E_K) 
          + I_app) / C_m #Unless we are doing IC, this has to stay this way
     @. dn = (Λ(v, V3, V4) * ((N∞(v, V3, V4) - n))) / τn
     @. dm = α_M(v, V7, V8, V9) * (1 - m) - β_M(v, V10, V11, V12) * m
     @. dh = α_H(v, V13, V14, V15) * (1 - h) - β_H(v, V16, V17, V18) * h
     @. dc = (C_0 + δ * (ICa(v, g_Ca, V1, V2, E_Ca)) - λ * c) / τc
     @. da = (α * c^a_n * (1 - a) - a) / τa #These were the old options
     @. db = (β * a^b_n * (1 - b) - b) / τb #These were the old options
     nothing
end
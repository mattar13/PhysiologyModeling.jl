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
     W = view(u, 11)

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
     dW = view(du, 11)

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
          V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18
     ) = extract_p0(p)

     @. dI_ext = I_app-I_ext
     @. dv = (ILeak(v, g_leak, E_leak) + 
          ICa(v, g_Ca, V1, V2, E_Ca) + IK(v, n, g_K, E_K) + ITREK(v, b, g_TREK, E_K) + INa(v, m, h, g_Na, E_Na) +
          IACh(v, e, g_ACh, k_ACh, E_ACh) + IGABA(v, i, g_GABA, k_GABA, E_Cl) + 
          I_app + W) / C_m #Unless we are doing IC, this has to stay this way
     @. dn = (Λ(v, V3, V4) * ((N∞(v, V3, V4) - n))) / τn
     @. dm = α_M(v, V7, V8, V9) * (1 - m) - β_M(v, V10, V11, V12) * m
     @. dh = α_H(v, V13, V14, V15) * (1 - h) - β_H(v, V16, V17, V18) * h
     @. dc = (C_0 + δ * (ICa(v, g_Ca, V1, V2, E_Ca)) - λ * c) / τc
     @. da = (α * c^a_n * (1 - a) - a) / τa #These were the old options
     @. db = (β * a^b_n * (1 - b) - b) / τb #These were the old options
     #@. da = (-α*(c^a_n)*a + (1-a))/τa     
     #@. db = (β * (1-a)^b_n * (1 - b) - b) / τb
     @. de = (ρe * Φe(v, VSe, V0e) - e) / τACh
     @. di = (ρi * Φi(v, VSi, V0i) - i) / τGABA
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
     W = view(u, 11)

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
     dW = view(du, 11)

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
          V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18
     ) = extract_p0(p)

     if stim_start < t < stim_stop
          stim_amp = I_app
     else
          stim_amp = 0.0
     end
     @. dI_ext = stim_amp - I_ext
     @. dv = (ILeak(v, g_leak, E_leak) + 
          ICa(v, g_Ca, V1, V2, E_Ca) + IK(v, n, g_K, E_K) + ITREK(v, b, g_TREK, E_K) + INa(v, m, h, g_Na, E_Na) +
          IACh(v, e, g_ACh, k_ACh, E_ACh) + IGABA(v, i, g_GABA, k_GABA, E_Cl) + 
          I_ext + W) / C_m
     @. dn = (Λ(v, V3, V4) * ((N∞(v, V3, V4) - n))) / τn
     @. dm = α_M(v, V7, V8, V9) * (1 - m) - β_M(v, V10, V11, V12) * m
     @. dh = α_H(v, V13, V14, V15) * (1 - h) - β_H(v, V16, V17, V18) * h
     @. dc = (C_0 + δ * (ICa(v, g_Ca, V1, V2, E_Ca)) - λ * c) / τc
     @. da = (α * c^a_n * (1 - a) - a) / τa #These were the old options
     @. db = (β * a^b_n * (1 - b) - b) / τb #These were the old options
     #@. da = (-α*(c^a_n)*a + (1-a))/τa     
     #@. db = (β * (1-a)^b_n * (1 - b) - b) / τb
     @. de = (ρe * Φe(v, VSe, V0e) - e) / τACh
     @. di = (ρi * Φi(v, VSi, V0i) - i) / τGABA
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
     W = view(u, 11)

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
     dW = view(du, 11)

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
          V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18
     ) = extract_p0(p)
     
     if isnothing(stim_start) && isnothing(stim_stop)
          @. dv = (VC-v)*k
     elseif stim_start < t < stim_stop
          @. dv = (VC-v)*k
     else
          @. dv = (hold-v)*k
     end
     @. dI_ext = (ILeak(v, g_leak, E_leak) + 
          ICa(v, g_Ca, V1, V2, E_Ca) + 
          IK(v, n, g_K, E_K) + ITREK(v, b, g_TREK, E_K) + INa(v, m, h, g_Na, E_Na) +
          IACh(v, e, g_ACh, k_ACh, E_ACh) + IGABA(v, i, g_GABA, k_GABA, E_Cl) + W) - I_ext
     @. dn = (Λ(v, V3, V4) * ((N∞(v, V3, V4) - n))) / τn
     @. dm = α_M(v, V7, V8, V9) * (1 - m) - β_M(v, V10, V11, V12) * m
     @. dh = α_H(v, V13, V14, V15) * (1 - h) - β_H(v, V16, V17, V18) * h
     @. dc = (C_0 + δ * (ICa(v, g_Ca, V1, V2, E_Ca)) - λ * c) / τc
     @. da = (α * c^a_n * (1 - a) - a) / τa #These were the old options
     @. db = (β * a^b_n * (1 - b) - b) / τb #These were the old options
     #@. da = (-α*(c^a_n)*a + (1-a))/τa     
     #@. db = (β * (1-a)^b_n * (1 - b) - b) / τb
     @. de = (ρe * Φe(v, VSe, V0e) - e) / τACh
     @. di = (ρi * Φi(v, VSi, V0i) - i) / τGABA
     @. dW = -W / τw

     nothing
end

function SAC_ODE_INH_EXC(du, u, p, t)
     I_ext = view(u, 1)
     v = view(u, 2)
     n = view(u, 3)
     m = view(u, 4)
     h = view(u, 5)
     c = view(u, 6)
     a = view(u, 7)
     b = view(u, 8)
     gE = view(u, 9)
     gI = view(u, 10)
     W = view(u, 11)

     dI_ext = view(du, 1)
     dv = view(du, 2)
     dn = view(du, 3)
     dm = view(du, 4)
     dh = view(du, 5)
     dc = view(du, 6)
     da = view(du, 7)
     db = view(du, 8)
     dgE = view(du, 9)
     dgI = view(du, 10)
     dW = view(du, 11)

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
          VSe, V0e, ρe,  g_ACh, k_ACh, E_ACh,  τACh, #Many of these will go unused
          VSi, V0i, ρi,  g_GABA, k_GABA, E_Cl, τGABA, #Many of these will go unused
          V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18
     ) = extract_p0(p)

     @. dI_ext = I_app-I_ext
     @. dv = (ILeak(v, g_leak, E_leak) + ICa(v, g_Ca, V1, V2, E_Ca) + IK(v, n, g_K, E_K) + ITREK(v, b, g_TREK, E_K) + INa(v, m, h, g_Na, E_Na) +
          + -gE*(v-E_ACh) #The conductance will be determined
          + -gI*(v-E_Cl) #The conductance will be determined
          + I_ext + W) / C_m
     @. dn = (Λ(v, V3, V4) * ((N∞(v, V3, V4) - n))) / τn
     @. dm = α_M(v, V7, V8, V9) * (1 - m) - β_M(v, V10, V11, V12) * m
     @. dh = α_H(v, V13, V14, V15) * (1 - h) - β_H(v, V16, V17, V18) * h
     @. dc = (C_0 + δ * (ICa(v, g_Ca, V1, V2, E_Ca)) - λ * c) / τc
     @. da = (α * c^a_n * (1 - a) - a) / τa #These were the old options
     @. db = (β * a^b_n * (1 - b) - b) / τb #These were the old options
     #@. da = (-α*(c^a_n)*a + (1-a))/τa     
     #@. db = (β * (1-a)^b_n * (1 - b) - b) / τb
     @. dgE = -gE
     @. dgI = -gI
     @. dW = -W / τw
     nothing
end

function SAC_ODE_Compartment(du, u, p, t; 
     gGAP = 0.01, 
     gK = [ 10.0, 8.125, 6.25, 4.375, 0.5], 
     #stim_starts = [2000.0, 1500.0, 1000.0, 500.0, 0.0], #Reverse the
     #stim_starts = [0.0, 500.0, 1000.0, 1500.0, 2000.0], #Forward list
     stim_starts = [0.0, 1000.0, 2000.0, 3000.0, 10000.0], #Reverse the
     stim_amp = 100.0,
     stim_dur = 500.0
)
     for i in axes(u, 1)
          dui = view(du, i, :)
          ui = view(u, i, :)
          p[7] = gK[i] #Change the parameter
          if i == 1
               SAC_ODE_STIM(dui, ui, p, t; stim_start = stim_starts[i], stim_stop = stim_starts[i]+stim_dur, stim_amp = stim_amp)
          else
               #println(500*(i-1))
               #Change the paramater gK
               SAC_ODE_STIM(dui, ui, p, t; stim_start = stim_starts[i], stim_stop = stim_starts[i]+stim_dur, stim_amp = stim_amp)
               dui[1] += gGAP * (u[i-1, 1] - u[i, 1]) #maybe make this more eff
          end
     end
end

#A more inline version
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
     W = view(u, :, 11)

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
     dW = view(du, :, 11)

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
          V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18
     ) = extract_p0(p)

     @. dI_ext = I_app
     @. dv = (ILeak(v, g_leak, E_leak) + 
          ICa(v, g_Ca, V1, V2, E_Ca) + IK(v, n, g_K, E_K) + ITREK(v, b, g_TREK, E_K) + INa(v, m, h, g_Na, E_Na) +
          IACh(v, e, g_ACh, k_ACh, E_ACh) + IGABA(v, i, g_GABA, k_GABA, E_Cl) + 
          I_ext + W) / C_m 
     @. dn = (Λ(v, V3, V4) * ((N∞(v, V3, V4) - n))) / τn
     @. dm = α_M(v, V7, V8, V9) * (1 - m) - β_M(v, V10, V11, V12) * m
     @. dh = α_H(v, V13, V14, V15) * (1 - h) - β_H(v, V16, V17, V18) * h
     @. dc = (C_0 + δ * (ICa(v, g_Ca, V1, V2, E_Ca)) - λ * c) / τc
     @. da = (α * c^a_n * (1 - a) - a) / τa #These were the old options
     @. db = (β * a^b_n * (1 - b) - b) / τb #These were the old options
     #@. da = (-α*(c^a_n)*a + (1-a))/τa     
     #@. db = (β * (1-a)^b_n * (1 - b) - b) / τb
     @. de = (ρe * Φe(v, VSe, V0e) - e) / τACh
     ∇α(de, e, MAP, t) #This takes alot of allocations. 
     @. di = (ρi * Φi(v, VSi, V0i) - i) / τGABA
     @. dW = -W / τw

     return
end

#Periodic stimulation of a x value
function SAC_PDE_STIM(du, u, MAP_p, t)
     #p[1] will be the cell map necessary for the equation to be run
     cell_map = MAP_p[1]
     #p[2] is the parameter set
     p = MAP_p[2]
     for i in axes(u, 1)
          if i == 61 && 5.0 < t < 100.0
               p[1] = 10.0
          else
               p[1] = 0.0
          end 
          dui = view(du, i, :)
          ui = view(u, i, :)
          SAC_ODE(dui, ui, p, t)
     end
     dE = view(du, :, 8)
     E = view(u, :, 8)
     ∇α(dE, E, cell_map, t)
end

DIFFUSION_NOISE(du, u, p, t) = du[:] .= 0.001

noise2D(du, u, p, t) = du[:, end] .= p[4]
function noise1D(du, u, p, t) 
     du[end] = p[4]
end

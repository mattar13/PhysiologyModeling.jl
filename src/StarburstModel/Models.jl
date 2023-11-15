function SAC_ODE(du, u, p, t)
     v = view(u, 1)
     n = view(u, 2)
     m = view(u, 3)
     h = view(u, 4)
     c = view(u, 5)
     a = view(u, 6)
     b = view(u, 7)
     e = view(u, 8)
     i = view(u, 9)
     W = view(u, 10)

     dv = view(du, 1)
     dn = view(du, 2)
     dm = view(du, 3)
     dh = view(du, 4)
     dc = view(du, 5)
     da = view(du, 6)
     db = view(du, 7)
     de = view(du, 8)
     di = view(du, 9)
     dW = view(du, 10)

     (I_app,
          C_m, g_W, τw, 
          g_leak, E_leak, 
          g_K, V3, V4, E_K, τn, 
          g_Ca, V1, V2,E_Ca, τc,
          g_Na, E_Na, 
          g_TREK,
          C_0, λ , δ,  
          α, τa, 
          β, τb, 
          VSe, ρe, V0e, g_ACh, k_ACh, E_ACh,  τACh,
          VSi, V0i, ρi,  g_GABA, k_GABA, E_Cl, τGABA,
          De, Di, 
          V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18
     ) = extract_p0(p)

     @. dv = (ILeak(v, g_leak, E_leak) + 
          ICa(v, g_Ca, V1, V2, E_Ca) + IK(v, n, g_K, E_K) + ITREK(v, b, g_TREK, E_K) + INa(v, m, h, g_Na, E_Na) +
          IACh(v, e, g_ACh, k_ACh, E_ACh) + IGABA(v, i, g_GABA, k_GABA, E_Cl) + 
          I_app + W) / C_m
     @. dn = (Λ(v, V3, V4) * ((N∞(v, V3, V4) - n))) / τn
     @. dm = α_M(v, V7, V8, V9) * (1 - m) - β_M(v, V10, V11, V12) * m
     @. dh = α_H(v, V13, V14, V15) * (1 - h) - β_H(v, V16, V17, V18) * h
     @. dc = (C_0 + δ * (ICa(v, g_Ca, V1, V2, E_Ca)) - λ * c) / τc
     @. da = (α * c^4 * (1 - a) - a) / τa
     @. db = (β * a^4 * (1 - b) - b) / τb
     @. de = (ρe * Φe(v, VSe, V0e) - e) / τACh
     @. di = (ρi * Φi(v, VSi, V0i) - i) / τGABA
     @. dW = -W / τw
     nothing
end

function SAC_ODE_NT_CLAMP(du, u, p, t)
     v = view(u, 1)
     n = view(u, 2)
     m = view(u, 3)
     h = view(u, 4)
     c = view(u, 5)
     a = view(u, 6)
     b = view(u, 7)
     e = view(u, 8)
     i = view(u, 9)
     W = view(u, 10)

     dv = view(du, 1)
     dn = view(du, 2)
     dm = view(du, 3)
     dh = view(du, 4)
     dc = view(du, 5)
     da = view(du, 6)
     db = view(du, 7)
     de = view(du, 8)
     di = view(du, 9)
     dW = view(du, 10)

     (I_app,
          C_m, g_W, τw, 
          g_leak, E_leak, 
          g_K, V3, V4, E_K, τn, 
          g_Ca, V1, V2,E_Ca, τc,
          g_Na, E_Na, 
          g_TREK,
          C_0, λ , δ,  
          α, τa, 
          β, τb, 
          VSe, ρe, V0e, g_ACh, k_ACh, E_ACh,  τACh,
          VSi, V0i, ρi,  g_GABA, k_GABA, E_Cl, τGABA,
          De, Di, 
          V7, V8, V9, V10, V11, V12, V13, V14, V15, V16, V17, V18
     ) = extract_p0(p)

     @. dv = (ILeak(v, g_leak, E_leak) + 
          ICa(v, g_Ca, V1, V2, E_Ca) + IK(v, n, g_K, E_K) + ITREK(v, b, g_TREK, E_K) + INa(v, m, h, g_Na, E_Na) +
          IACh(v, e, g_ACh, k_ACh, E_ACh) + IGABA(v, i, g_GABA, k_GABA, E_Cl) + 
          I_app + W) / C_m
     @. dn = (Λ(v, V3, V4) * ((N∞(v, V3, V4) - n))) / τn
     @. dm = α_M(v, V7, V8, V9) * (1 - m) - β_M(v, V10, V11, V12) * m
     @. dh = α_H(v, V13, V14, V15) * (1 - h) - β_H(v, V16, V17, V18) * h
     @. dc = (C_0 + δ * (ICa(v, g_Ca, V1, V2, E_Ca)) - λ * c) / τc
     @. da = (α * c^4 * (1 - a) - a) / τa
     @. db = (β * a^4 * (1 - b) - b) / τb
     @. de = ρe-e 
     @. di = ρi-i
     @. dW = -W / τw
     nothing
end

noise1D(du, u, p, t) = du[end] = 0.1

function ∇α(du, u, cell_map, t)
     flow_out = similar(cell_map.strength[:, 1]) #preallocate an array for flow_out
     for cellx in eachindex(u)
          connected_indices = cell_map.connected_indices[cellx]
          mul!(flow_out, -cell_map.strength[:, cellx], u[cellx])
          du[cellx] += sum(flow_out)
          for (i, celly) in enumerate(connected_indices)
               flow_in = cell_map.strength[celly, cellx] * u[cellx]
               du[celly] += flow_in
          end
          #we need to find out if this cell is a boundary cell and subtract the correct amount from it
          #if the cell is a single boundary we subtract one flow out
     end
end

function DIFFUSION_MODEL(du, u, p, t)
     du .= -u/540 #du decays over time
     if 500.0 < t < 2500.0
          #We want to add some diffusive material during a time range
          du[221] += 0.1
     end
     ∇α(du, u, p, t)#Diffusion occurs after
     #We should go through and decay the edges 
     
end

DIFFUSION_NOISE(du, u, p, t) = du[:] .= 0.001

function SAC_PDE(du, u, MAP_p, t)
     #p[1] will be the cell map necessary for the equation to be run
     cell_map = MAP_p[1]
     #p[2] is the parameter set
     p = MAP_p[300]
     for i in axes(u, 1)
          dui = view(du, i, :)
          ui = view(u, i, :)
          SAC_ODE(dui, ui, p, t)
     end
     dE = view(du, :, 8)
     E = view(u, :, 8)
     ∇α(dE, E, cell_map, t)
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

noise2D(du, u, p, t) = du[:, end] .= 0.1
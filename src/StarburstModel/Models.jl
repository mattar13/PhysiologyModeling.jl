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
     @. de = (ρe * Φe(v, VSe, V0e) - e) / τACh #ρe-e, #
     @. di = (ρi * Φi(v, VSi, V0i) - i) / τGABA #ρi-i, #
     @. dW = -W / τw
     nothing
end

noise1D(du, u, p, t) = du[end] = 0.1

function ∇α(du, u, p, t)
     cell_map = p
     du .= 0.0 #Set all changes to 0
     for cellx in eachindex(u)
          connected_indices = findall(x -> x != 0, cell_map.connections[:, cellx])
          flow_out = -(cell_map.strength[:, cellx].*cell_map.connections[:, cellx]) * u[cellx]
          du[cellx] += sum(flow_out)
          for (i, celly) in enumerate(connected_indices)
               flow_in = cell_map.strength[celly, cellx]*cell_map.connections[celly, cellx] * u[cellx]
               du[celly] += flow_in
          end
          #we need to find out if this cell is a boundary cell and subtract the correct amount from it
          #if the cell is a single boundary we subtract one flow out
     end
end

function SAC_eqs_PDE(du, u, p, t)
     


end
InhExc_p0_dict = Dict(
    "C_m"       => 80e-12,
    "g_leak"    => 4.2e-9,
    "E_leak"    => -55e-3,
    "g_Exc"     => 1.0, 
    "E_Exc"     => 0.0e-3, 
    "g_Inh"     => 1.0,
    "E_Inh"     => -60e-3,
)


InhExc_p0_keys = ["C_m", "g_leak", "E_leak", "g_Exc", "E_Exc", "g_Inh", "E_Inh"]
vals = extract_dict(InhExc_p0_dict, InhExc_p0_keys)
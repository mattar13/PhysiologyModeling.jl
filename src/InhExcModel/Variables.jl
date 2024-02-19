InhExc_p0_dict = Dict(
    "C_m"       => 0.08, #Scaleed at e-12
    "g_leak"    => 4.2, #Something about scaling makes me feel wreird e-9
    "E_leak"    => -55.0,
    "E_Exc"     => 0.0, 
    "E_Inh"     => -60.0,
)

InhExc_p0_keys = ["C_m", "g_leak", "E_leak", "E_Exc", "E_Inh"]
vals = extract_dict(InhExc_p0_dict, InhExc_p0_keys)
SAC_u0_dict = Dict(
     "v" =>  -64.23420980876107,   #1  v(t) = -63.6 
     "n" => 1.354100901552747e-5,  #2  n(t) = 0.000 
     "m" => 0.057910684157464414,  #3  m(t) = 0.062 
     "h" => 0.5691112964636721,    #4  h(t) = 0.550 
     "c" => 0.07729839999918353,   #5  c(t) = 0.085 
     "a" => 0.021826140913040696,  #6  a(t) = 0.026 
     "b" => 7.714363674074932e-6,  #7  b(t) = 0.000 
     "e" => 0.04675168624810658,   #8  e(t) = 0.000 
     "i" => 0.03895973860120606,   #9  i(t) = 0.000 
     "W" => 0.000                  #10 W(t) = 0.000 #Initial conditions
)

SAC_p0_dict = Dict(
    "I_app"     => 0.0,
    "C_m"       => 13.6,
    "g_W"       => 0.1,
    "τw"        => 800.0,
    "g_leak"    => 2.0,
    "E_leak"    => -70.0,
    "g_K"       => 4.0,
    "V3"        => -25.0,
    "V4"        => 7.0,
    "E_K"       => -90.0,
    "τn"        => 5.0,
    "g_Ca"      => 8.5,
    "V1"        => -20.0,
    "V2"        => 20.0,
    "E_Ca"      => 50.0,
    "τc"        => 2000.0,
    "g_Na"      => 2.0,
    "E_Na"      => 55.0,
    "g_TREK"    => 4.0,
    
    "C_0"       => 0.088,
    "λ"         => 2.702,
    "δ"         => 0.010503,
    "α"         => 625.0,
    "τa"        => 8300.0,
    "β"         => 34.0, #old value was 34.0
    "τb"        => 8300.0,
    "a_n"       => 4.0,
    "b_n"       => 4.0,

    "ρe"        => 6.0,
    "VSe"       => 0.2,
    "V0e"       => -40.0,
    "τACh"      => 540.0,
    "g_ACh"     => 0.215,
    "k_ACh"     => 0.1,
    "E_ACh"     => 0.0,
    "ρi"        => 5.0,
    "VSi"       => 0.2,
    "V0i"       => -40.0,
    "τGABA"     => 1000.0, #1000
    "g_GABA"    => 0.9,
    "k_GABA"    => 0.1,
    "E_Cl"      => -65.0,

    "De"        => 0.01,
    "Di"        => 0.01,
    "V7"        => 10.0,
    "V8"        => -40.0,
    "V9"        => 10.0,
    "V10"       => 4.0,
    "V11"       => -65.0,
    "V12"       => 18.0,
    "V13"       => 0.07,
    "V14"       => -65.0,
    "V15"       => 20.0,
    "V16"       => 1.0,
    "V17"       => -35.0,
    "V18"       => 10.0
)

keys_u0 = ["v", "n", "m", "h", "c", "a", "b", "e", "i", "W"]
GAP_keys_u0 = ["v", "n", "m", "h", "c", "a", "b", "W"] #This parameter set is for gap junctions only

vals_u0 = map(k -> SAC_u0_dict[k], keys_u0)
nt_u0 = NamedTuple{Symbol.(keys_u0) |> Tuple}(vals_u0)

keys_p0 = [
     "I_app",
     "C_m", "g_W", "τw", 
     "g_leak", "E_leak", 
     "g_K", "V3", "V4", "E_K", "τn", 
     "g_Ca", "V1", "V2","E_Ca", "τc",
     "g_Na", "E_Na", 
     "g_TREK",
     
     "C_0", "λ" , "δ",  
     "α", "τa", 
     "β", "τb", 
     "a_n", "b_n",

     "VSe", "ρe", "V0e", "g_ACh", "k_ACh", "E_ACh",  "τACh",
     "VSi", "V0i", "ρi",  "g_GABA", "k_GABA", "E_Cl", "τGABA",

     "De", "Di", 
     "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18"
]

vals_p0 = map(k -> SAC_p0_dict[k], keys_p0) 

nt_p0 = NamedTuple{Symbol.(keys_p0) |> Tuple}(vals_p0)

extract_dict(d::Dict{String, Float64}, keys) = map(k -> d[k], keys)
extract_p0(d::Dict{String, Float64})  = extract_dict(d, keys_p0)
extract_p0(d::NamedTuple) = d
extract_p0(d) = d

extract_u0(d::Dict{String, Float64}) = extract_dict(d, keys_u0)
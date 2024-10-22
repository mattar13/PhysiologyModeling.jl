extract_dict(d::Dict{String, Float64}, keys) = map(k -> d[k], keys)
extract_dict(d::Dict{String, Vector{Float64}}, keys) = hcat(map(k -> d[k], keys)...)

function SAC_u0_dict(;mode = :ODE, n_cells = 1) 
     u0_dict = Dict(
          "I_ext" => 0.0,               #1  I_ext(t) = 0.00 This is the added current
          "v" =>  -64.23420980876107,   #2  v(t) = -63.6 
          "n" => 1.354100901552747e-5,  #3  n(t) = 0.000 
          "m" => 0.057910684157464414,  #4  m(t) = 0.062 
          "h" => 0.5691112964636721,    #5  h(t) = 0.550 
          "c" => 0.07729839999918353,   #6  c(t) = 0.085 
          "a" => 0.9768639635812157,  #7  a(t) = 0.026 
          "b" => 0.0003635613789261572,  #8  b(t) = 0.000 
          "e" => 0.04675168624810658,   #9  e(t) = 0.000 
          "i" => 0.03895973860120606,   #10 i(t) = 0.000 
          "g" => 0.0,
          "d" => 0.0, 
          "q" => 0.0,
          "W" => 0.000                  #11 W(t) = 0.000 #Initial conditions
     )
     if n_cells > 1
          u0_dict_vector = Dict{String, Vector{Float64}}()
          for (k,v) in u0_dict
               u0_dict_vector[k] = fill(v, n_cells)
          end
          return u0_dict_vector
     elseif mode == :ODE
          return u0_dict
     elseif mode == :GLUT
          u0_dict["g"] = 0.0
          u0_dict["q"] = 0.0 #Modulatory g-protein deactivating gCa
     end
end
#
function SAC_p0_dict(;keyset = :DEFAULT)
     base_dict = Dict{String, Union{Real, Vector}}(
          "I_app"     => 0.0,
          "VC"        => 0.0, #For voltage clamp
          "gGAP"      => 0.2, #for gap junction conductance
          "C_m"       => 13.6,
          "g_W"       => 0.075,
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

          "k_SYT"     => 1.0e-8,
          "n_SYT"     => 10.0,

          "ρe"        => 1.0,#6.0,
          "VSe"       => 0.2,
          "V0e"       => -40.0,
          "τACh"      => 124.0,
          "g_ACh"     => 2.15,
          "k_ACh"     => 0.1, #Half maximal concentration is 100 uM
          "E_ACh"     => 0.0,
          
          "ρi"        => 1.0,#5.0,
          "VSi"       => 0.2,
          "V0i"       => -40.0,
          "τGABA"     => 100.0, #1000
          "g_GABA"    => 0.9,
          "k_GABA"    => 0.1,
          "E_Cl"      => -65.0,

          #These parameters were estimated from Smith and Howe et al. 
          "g_GLUT"    => 0.15, #Single channel conductance is 4-15pS. Assuming 100 channels per cell. value is in nS
          "k_GLUT"    => 0.007, #halfmax activation occurs at 7 uM
          "E_GLUT"    => 0.0, #NMDA receptors are non-selective
          
          "τd"        => 5000.0, #decay of dopamine
          "g_n"       => 2.0, #number of glutamate moleculaes needed to activte the g protein response
          "γg"        => 6.18, #rate of Gq activation by G protein 
          "τq"        => 500.0, #time needed for Gq to return to baseline

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
          "V18"       => 10.0, 
          "stim_start" => -Inf, #Set it to this to always be on
          "stim_stop" => Inf
     )
     if keyset == :DEFAULT
          return base_dict
     elseif keyset == :dynamical_analysis
          base_dict["g_TREK"] = 0.0
          base_dict["g_ACh"] = 0.0
          base_dict["g_GABA"] = 0.0
          base_dict["g_GLUT"] = 0.0
          return base_dict
     end
end

function extract_p0(d::Dict{String, Union{Real, Vector}}; mode = :ODE)
     keys_p0 = [
          "I_app", "VC", "gGAP",
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

          "k_SYT", "n_SYT",

          "ρe",  "g_ACh", "k_ACh", "E_ACh",  "τACh",
          "ρi",  "g_GABA", "k_GABA", "E_Cl", "τGABA",

          "g_GLUT", "k_GLUT", "E_GLUT", 
          
          "τd",
          "γg", "g_n", "τq",

          "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18",
          "stim_start", "stim_stop"
     ]  
     return map(k -> d[k], keys_p0)
end

function extract_u0(d::Dict{String, T}; mode = :ODE) where T
     if mode == :ODE || mode == :PDE
          keys_u0 = ["I_ext", "v", "n", "m", "h", "c", "a", "b", "e", "i", "g", "d", "q", "W"]
     elseif mode == :KEYS

     elseif mode == :DynamicalAnalysis
          keys_u0 = ["v", "n", "m", "h", "c", "a", "b"] #This parameter set is for gap junctions only
     end
     extract_dict(d, keys_u0)
end

function par_idx(par::String)
     keys_p0 =keys_p0 = [
          "I_app", "VC", "gGAP",
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

          "k_SYT", "n_SYT",

          "ρe",  "g_ACh", "k_ACh", "E_ACh",  "τACh",
          "ρi",  "g_GABA", "k_GABA", "E_Cl", "τGABA",

          "g_GLUT", "k_GLUT", "E_GLUT", 
          
          "τd",
          "γg", "g_n", "τq",

          "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18",
          "stim_start", "stim_stop"
     ]  
     findfirst(keys_p0 .== par)
end
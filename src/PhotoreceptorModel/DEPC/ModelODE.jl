using Revise
using PhysiologyModeling
using PyPlot; PyPlot.pygui(true)

#%% Declare the variables ════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════#
@parameters t #One time dimension
dt = 0.001
tspan = (0.0, 10.0) #Define the time domain
tstops = tspan[1]:dt:tspan[2]

Dt = Differential(t) #First time derivative
#Voltage properties of the model
@variables hv(t) I_C(t) V(t) mKV(t) hKV(t) mCa(t) mKCa(t)

#Rxn network variables and parameters
@variables R(t) Ra(t) 
@variables Trαβγ(t) Trα(t) Trαi(t) Trβγ(t) 
@variables PDE(t) PDEa(t) 
@variables Ca(t) CaB(t) cGMP(t) GMP(t)

#Calcium Network variables and parameters
@variables _Ca_s(t) _Ca_f(t) _CaB_ls(t) _CaB_hs(t) _CaB_lf(t) _CaB_hf(t)

@variables C1(t) C2(t) O1(t) O2(t) O3(t) #This is a matrix

@variables GL(t)

@variables I_ALL(t)

@parameters gH EH
@parameters k_Ca τGL
@parameters cm gDARK gL EL gKV EK gCa ECa gCl ECl gKCa
@parameters gH EH
@parameters αC K1r K2f K2r K3f K3r K4 K5 K6f K6r
@parameters AMAX Kc b γCa C0 er
@parameters DCa Lb1 Lb2 Hb1 Hb2
@parameters Bl Bh
@parameters Cae J_ex J_ex2 K_ex K_ex2
@parameters F V1 S1 δ Lb1 Lb2 Hb1 Hb2 Bl Bh DCa _Ca_O #Calcium network parameters
@parameters F V1 V2 S1 δ
@parameters α C_m g_L E_L g_KV E_K g_Ca E_Ca g_Cl E_Cl g_KCa J_hv #Current parameters
@parameters J_ex J_ex2 Cae K_ex K_ex2 
@parameters μ ω

# Make the auxillary equations ════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════#
#Current equations
I_PHOTO(v, cGMP) = -gDARK * J(cGMP) * (1.0 - exp((v - 8.5) / 17.0))
I_L(v) = gL * (v - EL)
I_KV(v) = gKV * mKV^3 * hKV * (v - EK)
I_Ca(v) = gCa * mCa^4 * hCa(v) * (v - ECa)
I_Cl(v, _Ca_s) = gCl * mCl(_Ca_s) * (v - ECl)
I_KCa(v, _Ca_s) = gKCa * mKCa^2 * mKCas(_Ca_s) * (v - EK)
I_ex(v, _Ca_s) = J_ex * exp(-(v - 14) / 70) * (_Ca_s - Cae) / (_Ca_s - Cae + K_ex)
I_ex2(_Ca_s) = J_ex2 * (_Ca_s - Cae) / (_Ca_s - Cae + K_ex2)
I_h(v) = gH * (O1 + O1 + O3) * (v - EH)

J(cGMP) = (cGMP^3)/(cGMP^3 + 10^3)
A(Ca) = AMAX/(1.0 + (Ca / Kc)^4)

#KV equations
αmKV(v) = (5 * (100 - v)) / (exp((100 - v) / 42) - 1)
βmKV(v) = 9 * exp(-(v - 20) / 40)
αhKV(v) = 0.15 * exp(-v / 22)
βhKV(v) = 0.4125 / (exp((10 - v) / 7) + 1)

#Calcium equations
αmCa(v) = (3 * (80 - v)) / (exp((80 - v) / 25) - 1)
βmCa(v) = 10 / (1 + exp((v + 38) / 7))
hCa(v) = exp((40 - v) / 18) / (1 + exp((40 - v) / 18))

αmKCa(v) = (15 * (80 - v)) / (exp((80 - v) / 40) - 1)
βmKCa(v) = 20 * exp(-v / 35)
mKCas(_Ca_s) = _Ca_s / (_Ca_s + 0.3)

#Calcium activated chloride
mCl(_Ca_s) = 1 / (1 + exp((0.37 - _Ca_s) / 0.9))

#Calcium activated potassium
αmKCa(v) = (15 * (80 - v)) / (exp((80 - v) / 40) - 1)
βmKCa(v) = 20 * exp(-v / 35)
mKCas(_Ca_s) = _Ca_s / (_Ca_s + 0.3)

#Hyperpolarising current
αh(v) = 8 / (exp((v + 78) / 14) + 1)
βh(v) = 18 / (exp(-(v + 9) / 19) + 1)

#Hyperpolarising current
αh(v) = 8 / (exp((v + 78) / 14) + 1)
βh(v) = 18 / (exp(-(v + 9) / 19) + 1)

hT(v) = [
     -4*αh(v) βh(v) 0.0 0.0 0.0
     4*αh(v) -(3*αh(v) + βh(v)) 2*βh(v) 0.0 0.0
     0.0 3*αh(v) -(2*αh(v)+2*βh(v)) 3*βh(v) 0.0
     0.0 0.0 2*αh(v) -(αh(v) + 3*βh(v)) 4*βh(v)
     0.0 0.0 0.0 αh(v) -4*βh(v)
]

fhT(v, C1, C2, O1, O2, O3) = hT(v) * [C1, C2, O1, O2, O3]

# Define the equations ════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════#
@named PC_ALG = ODESystem([
     Dt(hv) ~ 0.0,
     #Set up the voltage equation
	Dt(V) ~ -(
          I_L(V) + I_KV(V) + I_Ca(V) + 
          I_Cl(V, _Ca_s) + I_KCa(V, _Ca_s) + 
          I_PHOTO(V, cGMP) + 
          I_ex(V, _Ca_s) + I_ex2(_Ca_s) + 
          I_h(V)
     )/cm,
     I_ALL ~ (I_L(V) + I_KV(V) + I_Ca(V) + 
          I_Cl(V, _Ca_s) + I_KCa(V, _Ca_s) + 
          I_PHOTO(V, cGMP) + 
          I_ex(V, _Ca_s) + I_ex2(_Ca_s) + 
          I_h(V)),
     #This is the reaction network
	Dt(R) ~ K1r*Ra - αC*R*hv,
	Dt(Ra) ~ αC*R*hv - K1r*Ra,
	Dt(Trαβγ) ~ K2r*Trαi*Trβγ - K2f*Ra*Trαβγ,
	Dt(Trα) ~ K2f*Ra*Trαβγ - K3f*PDE*Trα,
	Dt(Trβγ) ~ K2f*Ra*Trαβγ - K2r*Trαi*Trβγ,
	Dt(PDE) ~ K3r*PDEa - K3f*PDE*Trα,
	Dt(PDEa) ~ K3f*PDE*Trα - K3r*PDEa,
	Dt(Trαi) ~ K3r*PDEa - K2r*Trαi*Trβγ,
	Dt(cGMP) ~ (AMAX*GMP) / (1.0 + (Ca / Kc)^4) - K4*PDEa*cGMP,
	Dt(GMP) ~ (-AMAX*GMP) / (1.0 + (Ca / Kc)^4) + K4*PDEa*cGMP,
	Dt(Ca) ~ K6r*CaB + (b*gDARK*(cGMP^3)) / (1000 + cGMP^3) - K6f*(er - CaB)*Ca - γCa*(1 + (-C0) / Ca)*Ca,
	Dt(CaB) ~ K6f*(er - CaB)*Ca - K6r*CaB,
     #Intracellular calcium system
     Dt(_Ca_s) ~ -((I_Ca(V) + I_ex(V, _Ca_s) + I_ex2(_Ca_s)) / (2 * F * V1)) * 10e-6 - DCa * (S1 / (δ * V1)) * (_Ca_s - _Ca_f) - Lb1 * _Ca_s * (Bl - _CaB_ls) + Lb2 * _CaB_ls - Hb1 * _Ca_s * (Bh - _CaB_hs) + Hb2 * _CaB_hs,
     Dt(_Ca_f) ~ DCa * (S1 / (δ * V1)) * (_Ca_s - _Ca_f) - Lb1 * _Ca_f * (Bl - _CaB_lf) + Lb2 * _CaB_lf - Hb1 * _Ca_f * (Bh - _CaB_hf) + Hb2 * _CaB_hf,
     Dt(_CaB_ls) ~ Lb1 * _Ca_s * (Bl - _CaB_ls) - Lb2 * _CaB_ls,
     Dt(_CaB_hs) ~ Hb1 * _Ca_s * (Bh - _CaB_hs) - Hb2 * _CaB_hs,
     Dt(_CaB_lf) ~ Lb1 * _Ca_f * (Bl - _CaB_lf) - Lb2 * _CaB_lf,
     Dt(_CaB_hf) ~ Hb1 * _Ca_f * (Bh - _CaB_hf) - Hb2 * _CaB_hf,
     #These are non-linear voltage properties of the kv and ca channels
     Dt(mKV) ~ αmKV(V) * (1 - mKV) - βmKV(V) * mKV,
     Dt(hKV) ~ αhKV(V) * (1 - hKV) - βhKV(V) * hKV,
     Dt(mCa) ~ αmCa(V) * (1 - mCa) - βmCa(V) * mCa,
     #This is the non-linear term for the Calcium activate potassium
     Dt(mKCa) ~ αmKCa(V) * (1 - mKCa) - βmKCa(V) * mKCa,
     #Hyperpolarizing activated current
     Dt(C1) ~ fhT(V, C1, C2, O1, O2, O3)[1], #This transforms the hMAT variable
     Dt(C2) ~ fhT(V, C1, C2, O1, O2, O3)[2], #This transforms the hMAT variable
     Dt(O1) ~ fhT(V, C1, C2, O1, O2, O3)[3], #This transforms the hMAT variable
     Dt(O2) ~ fhT(V, C1, C2, O1, O2, O3)[4], #This transforms the hMAT variable
     Dt(O3) ~ fhT(V, C1, C2, O1, O2, O3)[5], #This transforms the hMAT variable
     #Glutamate release
     #Dt(GL) ~ (-k_Ca * I_Ca(V) - GL) / τGL
])

# Set the initial conditions and parameters ════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════#
u0 = [
	hv => 0.0      #Photons/μm
	V => -35.0     #Membrane Voltage
	mKV => 0.430   #Rectifying potassium variable 1
	hKV => 0.999   #Rectifying potassium variable 2
	mCa => 0.436   #Voltage gated calcium current
	mKCa => 0.642  #Calcium activated potassium channel
	
	#Hyperpolarising currents
	C1 => 0.646
     C2 => 0.298
     O1 => 0.0517
     O2 => 0.00398
    	O3 => 0.000115

	#Phototransduction cascade
	R => 10000.0
	Ra => 0.0
	Trαβγ => 1000.0
	Trα => 0.0
	Trαi => 0.0
	Trβγ => 0.0
	PDE => 100.0
	PDEa => 0.0
	cGMP => 2.0
	GMP => 0.0

	Ca => 0.3 #Initial calcium
	CaB => 272.72727272727286 #Calcium Buffers

	#Calcium binding proteins
	_Ca_s => 0.0966 #Free calcium
 	_Ca_f => 0.0966
 	_CaB_ls => 80.929
 	_CaB_hs => 29.068
	_CaB_lf => 80.929
	_CaB_hf => 29.068
	#GL => 0.0

     #Currents for viewing 
     I_ALL => 0.0
];
p0 = [
	αC => 0.6016      #Collecting area    DEFAULT: 0.6016
	K1r => 50.0 #Regeneration of Rhodopsin
	K2f => 0.5 #Transducin activation by Rhodopsin
	K2r => 2.5 #Regeneration of Transducin
	K3f => 0.2
	K3r => 5.0
	K4 => 1.4
	K5 => 0.8
	gDARK => 5040.0
	cm => 0.02
	gL => 0.35   #Leaky Conductance       DEFAULT:  0.3335 nS
	EL => -77.0 
	gKV => 2.0   #KV conductance          DEFAULT:  2.0 nS 
	EK => -74.0  #Potassium Reversel      DEFAULT: -74.0 mV
	AMAX => 65.6
	b => 0.25
	γCa => 50.0
	C0 => 0.1      #Prevent calcium from falling below this limit
	er => 500      #This is the limits of calcium buffering DEFAULT 500uM
 	K6r => 0.2     #Buffering of calcium    DEFAULT: 0.2 s^1 μM^-1 
 	K6f => 0.8     #Unbuffering of calcium  DEFAULT: 0.8 s^-1
	Kc => 0.1      #Dimensionless reverse   DEFAULT: 0.1
	DCa => 6e-8    #                        DEFAULT: 6.0e-8 dm^2 
	Lb1 => 0.4     #                        DEFAULT: 0.4 s^-1
	Lb2 => 0.2     #                        DEFAULT: 0.2 s^-1
	Hb1 => 100     #                        DEFAULT: 100 s^-1 μM^-1
	Hb2 => 90      #                        DEFAULT: 90 s^-1
	Bl => 500      #                        DEFAULT: 500 μM
	Bh => 300      #                        DEFAULT: 300 μm
	gCa => 0.7     #Calicum maximal cond    DEFAULT: 0.7 nS
	ECa => -12.5
	Cae => 0.01    #                        DEFAULT: 0.01 μM
	J_ex => 20.0   #                        DEFAULT: 20.0 pA
	J_ex2 => 20.0  #                        DEFAULT: 20.0 pA
	K_ex => 2.3    #                        DEFAULT: 2.3 μM 
	K_ex2 => 0.5   #                        DEFAULT: 0.5 μM
	S1 => 3.142e-8 #const                   DEFAULT: 3.142e-8 dm^2
	δ => 3e-5      #Calcium area            DEFAULT: 3.0e-5
	F => 9.648e4   #This is Faraday const   DEFAULT: 9.648e4 C mol^-1
	V1 => 3.812e-13 #A second const		DEFAULT: 3.812e-13 dm^3
	V2 => 5.236e-13 #A third const		DEFAULT: 5.236e-13 dm^3
     gCl => 2.0   #Chloride channel cond.  DEFAULT: 2.0 nS
     ECl => -20   #Cl Reversal             DEFAULT: -20.0 mV 
     gKCa => 5.0  #Calcium activated k     DEFAULT: 5.0 nS
     gH => 3.0     #h channels conductance  DEFAULT: 3.0 nS
     EH => -32.0   #h channel reversal      DEFAULT: -32.0 mV
];
#Make the model
PC_ODE = structural_simplify(PC_ALG)
stable_prob = ODEProblem(PC_ODE, u0, tspan, p0) 
#make a stable solution
stable_u = solve(stable_prob, save_start = false, save_everystep = false, save_end = true)
prob = remake(stable_prob, u0 = stable_u.u[end])

# Make the light stimulus ════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════#
stim_start = 2.0
stim_end = 2.01
stimulus = zeros(length(tstops))
stimulus[stim_start .< tstops .< stim_end] .= 1000

affect!(integrator) = integrator.u[1] = stimulus[floor(Int64, (integrator.t / dt) + 1)]
cb = PresetTimeCallback(tstops, affect!)

#%% Solve the model ════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════#
sol = solve(prob, Rodas5(), callback = cb, tstops = tstops)

#%% Plot the channel conductances ════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════#
Time = sol.t #Extract the time variables

fig3, axs = plt.subplots(3,1)
hvt = sol.(Time, idxs = hv)
axs[1].plot(Time, hvt)

Vt = sol.(Time, idxs = V)
cGMPt = sol.(Time, idxs = cGMP)
I_ALLt = sol.(Time, idxs = I_ALL)
#axs[2].plot(Time, cGMPt, label = "cGMP")
axs[2].plot(Time, I_ALLt, label = "Dark Current")
axs[3].plot(Time, Vt, label = "Voltage")

#%% Plot the Calcium buffer system  ════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════#
fig1, axs = plt.subplots(5,1)
axs[1].plot(Time, hvt)
Cat = sol.(Time, idxs = Ca)
CaBt = sol.(Time, idxs = CaB)
axs[2].plot(Time, Cat, label = "Ca")
axs[2].plot(Time, CaBt, label = "CaB")

_Ca_st = sol.(Time, idxs = _Ca_s)
_Ca_ft = sol.(Time, idxs = _Ca_f)
axs[3].plot(Time, _Ca_st, label = "[Ca]s")
axs[3].plot(Time, _Ca_ft, label = "[Ca]f")
axs[3].legend()

_CaB_lst = sol.(Time, idxs = _CaB_ls)
_CaB_hst = sol.(Time, idxs = _CaB_hs)
axs[4].plot(Time, _Ca_st, label = "[Ca]s")
axs[4].plot(Time, _CaB_lst, label = "[CaB]ls")
axs[4].plot(Time, _CaB_hst, label = "[CaB]hs")
axs[4].legend()

_CaB_lft = sol.(Time, idxs = _CaB_lf)
_CaB_hft = sol.(Time, idxs = _CaB_hf)
axs[5].plot(Time, _Ca_ft, label = "[Ca]f")
axs[5].plot(Time, _CaB_lft, label = "[CaB]lf")
axs[5].plot(Time, _CaB_hft, label = "[CaB]hf")
axs[5].legend()

#%% Plot the Phototransduction cascade  ════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════#
fig, axs = plt.subplots(5,1)
axs[1].plot(Time, hvt)
#axs[1].plot(Time, It)
axs[1].set_ylabel("Phot")
axs[1].xaxis.set_visible(false)

Rt = sol.(Time, idxs = R)
Rat = sol.(Time, idxs = Ra)
axs[2].plot(Time, Rt, label = "R")
axs[2].plot(Time, Rat, label = "R*")
axs[2].set_ylabel("Rh")
axs[2].xaxis.set_visible(false)
axs[2].legend()

Trαβγt = sol.(Time, idxs = Trαβγ)
Trαt = sol.(Time, idxs = Trα)
Trαit = sol.(Time, idxs = Trαi)
Trβγt = sol.(Time, idxs = Trβγ)
axs[3].plot(Time, Trαβγt, label = "Gnαβγ")
axs[3].plot(Time, Trαt, label = "Gnαt*")
axs[3].plot(Time, Trαit, label = "Gnαt")
axs[3].plot(Time, Trβγt, label = "Gnβγ")
axs[3].set_ylabel("GNAT")
axs[3].legend()
axs[3].xaxis.set_visible(false)

PDEt = sol.(Time, idxs = PDE)
PDEat = sol.(Time, idxs = PDEa)
axs[4].plot(Time, PDEt)
axs[4].plot(Time, PDEat)
axs[4].set_ylabel("PDE")
axs[4].xaxis.set_visible(false)

cGMPt = sol.(Time, idxs = cGMP)
GMPt = sol.(Time, idxs = GMP)
axs[5].plot(Time, cGMPt, label = "cGMP")
axs[5].plot(Time, GMPt, label = "GMP")
axs[5].set_ylabel("cGMP")
axs[5].legend()

#%% lets make an intensity response function from the model
fig4, axs = plt.subplots(3,1)
for Φ in 10 .^ (-4.0:1.0:5.0)
     println("Stimulating at light intensity: $Φ")
     # Make the light stimulus ════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════#
     stim_start = 2.0
     stim_end = 2.01
     stimulus = zeros(length(tstops))
     stimulus[stim_start .< tstops .< stim_end] .= Φ

     affect!(integrator) = integrator.u[1] = stimulus[floor(Int64, (integrator.t / dt) + 1)]
     cb = PresetTimeCallback(tstops, affect!)

     #%% Solve the model ════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════#
     sol_i = solve(prob, Rodas5(), callback = cb, tstops = tstops)
     Time = sol_i.t #Extract the time variables

     hvt = sol_i.(Time, idxs = hv)
     axs[1].plot(Time, hvt)

     I_ALLt = sol_i.(Time, idxs = I_ALL)
     axs[2].plot(Time, I_ALLt, label = "Dark Current")
     
     Vt = sol_i.(Time, idxs = V)
     axs[3].plot(Time, Vt, label = "Voltage")
end
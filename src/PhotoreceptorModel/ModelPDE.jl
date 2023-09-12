using PyPlot; PyPlot.pygui(true)
using DifferentialEquations, ModelingToolkit
using MethodOfLines, DomainSets
using Symbolics: scalarize
using Latexify

#%% Declare the variables ════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════#
@parameters x y t #One time dimension
Dt = Differential(t) #First time derivative
Dx = Differential(x)
Dy = Differential(y)
Dxx = Differential(x)^2
Dyy = Differential(y)^2

dt = 0.01; dx = 0.1; dy = 0.1
xmin = ymin = tmin = 0.0
tmax = 10.0
xmax = ymax = 1.0
tspan = (tmin, tmax) #Define the time domain
tstops = tspan[1]:dt:tspan[2]

#Voltage properties of the model
begin
	@variables hv(..) I_C(..) V(..) mKV(..) hKV(..) mCa(..) mKCa(..)

	#Rxn network variables and parameters
	@variables R(..) Ra(..) 
	@variables Trαβγ(..) Trα(..) Trαi(..) Trβγ(..) 
	@variables PDE(..) PDEa(..) 
	@variables Ca(..) CaB(..) cGMP(..) GMP(..)


	#Calcium Network variables and parameters
	@variables _Ca_s(..) _Ca_f(..) _CaB_ls(..) _CaB_hs(..) _CaB_lf(..) _CaB_hf(..)

	@variables C1(..) C2(..) O1(..) O2(..) O3(..) #This is a matrix

	@variables GL(..)

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
	@parameters gGAP
end
# Make the auxillary equations ════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════#
begin
	#Current equations
	I_PHOTO(v, cGMP) = -gDARK * J(cGMP) * (1.0 - exp((v - 8.5) / 17.0))
	I_L(v) = gL * (v - EL)
	I_KV(v) = gKV * mKV(x,y,t)^3 * hKV(x,y,t) * (v - EK)
	I_Ca(v) = gCa * mCa(x,y,t)^4 * hCa(v) * (v - ECa)
	I_Cl(v, _Ca_s) = gCl * mCl(_Ca_s) * (v - ECl)
	I_KCa(v, _Ca_s) = gKCa * mKCa(x,y,t)^2 * mKCas(_Ca_s) * (v - EK)
	I_ex(v, _Ca_s) = J_ex * exp(-(v - 14) / 70) * (_Ca_s - Cae) / (_Ca_s - Cae + K_ex)
	I_ex2(_Ca_s) = J_ex2 * (_Ca_s - Cae) / (_Ca_s - Cae + K_ex2)
	I_h(v) = gH * (O1(x,y,t) + O1(x,y,t) + O3(x,y,t)) * (v - EH)

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
	∇²(v) = Dxx(v) + Dyy(v) #Voltage crossing the gap junction
end
# Define the equations ════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════#
PDEeqs = [
	Dt(hv(x,y,t)) ~ 0.0,
	#Set up the voltage equation
	Dt(V(x,y,t)) ~ -(
		-gGAP * ∇²(V(x, y, t)) +
		I_L(V(x,y,t)) + I_KV(V(x,y,t)) + I_Ca(V(x,y,t)) + 
		I_Cl(V(x,y,t), _Ca_s(x,y,t)) + I_KCa(V(x,y,t), _Ca_s(x,y,t)) + 
		I_PHOTO(V(x,y,t), cGMP(x,y,t)) + 
		I_ex(V(x,y,t), _Ca_s(x,y,t)) + I_ex2(_Ca_s(x,y,t)) + 
		I_h(V(x,y,t))
	)/cm,
	#This is the reaction network
	Dt(R(x,y,t)) ~ K1r*Ra(x,y,t) - αC*R(x,y,t)*hv(x,y,t),
	Dt(Ra(x,y,t)) ~ αC*R(x,y,t)*hv(x,y,t) - K1r*Ra(x,y,t),
	Dt(Trαβγ(x,y,t)) ~ K2r*Trαi(x,y,t)*Trβγ(x,y,t) - K2f*Ra(x,y,t)*Trαβγ(x,y,t),
	Dt(Trα(x,y,t)) ~ K2f*Ra(x,y,t)*Trαβγ(x,y,t) - K3f*PDE(x,y,t)*Trα(x,y,t),
	Dt(Trβγ(x,y,t)) ~ K2f*Ra(x,y,t)*Trαβγ(x,y,t) - K2r*Trαi(x,y,t)*Trβγ(x,y,t),
	Dt(PDE(x,y,t)) ~ K3r*PDEa(x,y,t) - K3f*PDE(x,y,t)*Trα(x,y,t),
	Dt(PDEa(x,y,t)) ~ K3f*PDE(x,y,t)*Trα(x,y,t) - K3r*PDEa(x,y,t),
	Dt(Trαi(x,y,t)) ~ K3r*PDEa(x,y,t) - K2r*Trαi(x,y,t)*Trβγ(x,y,t),
	Dt(cGMP(x,y,t)) ~ (AMAX*GMP(x,y,t)) / (1.0 + (Ca(x,y,t) / Kc)^4) - K4*PDEa(x,y,t)*cGMP(x,y,t),
	Dt(GMP(x,y,t)) ~ (-AMAX*GMP(x,y,t)) / (1.0 + (Ca(x,y,t) / Kc)^4) + K4*PDEa(x,y,t)*cGMP(x,y,t),
	Dt(Ca(x,y,t)) ~ K6r*CaB(x,y,t) + (b*gDARK*(cGMP(x,y,t)^3)) / (1000 + cGMP(x,y,t)^3) - K6f*(er - CaB(x,y,t))*Ca(x,y,t) - γCa*(1 + (-C0) / Ca(x,y,t))*Ca(x,y,t),
	Dt(CaB(x,y,t)) ~ K6f*(er - CaB(x,y,t))*Ca(x,y,t) - K6r*CaB(x,y,t),
	#Intracellular calcium system
	Dt(_Ca_s(x,y,t)) ~ -((I_Ca(V(x,y,t)) + I_ex(V(x,y,t), _Ca_s(x,y,t)) + I_ex2(_Ca_s(x,y,t))) / (2 * F * V1)) * 10e-6 - DCa * (S1 / (δ * V1)) * (_Ca_s(x,y,t) - _Ca_f(x,y,t)) - Lb1 * _Ca_s(x,y,t) * (Bl - _CaB_ls(x,y,t)) + Lb2 * _CaB_ls(x,y,t) - Hb1 * _Ca_s(x,y,t) * (Bh - _CaB_ls(x,y,t)) + Hb2 * _CaB_hs(x,y,t),
	Dt(_Ca_f(x,y,t)) ~ DCa * (S1 / (δ * V1)) * (_Ca_s(x,y,t) - _Ca_f(x,y,t)) - Lb1 * _Ca_f(x,y,t) * (Bl - _CaB_lf(x,y,t)) + Lb2 * _CaB_lf(x,y,t) - Hb1 * _Ca_f(x,y,t) * (Bh - _CaB_hf(x,y,t)) + Hb2 * _CaB_hf(x,y,t),
	Dt(_CaB_ls(x,y,t)) ~ Lb1 * _Ca_s(x,y,t) * (Bl - _CaB_ls(x,y,t)) - Lb2 * _CaB_ls(x,y,t),
	Dt(_CaB_hs(x,y,t)) ~ Hb1 * _Ca_s(x,y,t) * (Bh - _CaB_hs(x,y,t)) - Hb2 * _CaB_hs(x,y,t),
	Dt(_CaB_lf(x,y,t)) ~ Lb1 * _Ca_f(x,y,t) * (Bl - _CaB_lf(x,y,t)) - Lb2 * _CaB_lf(x,y,t),
	Dt(_CaB_hf(x,y,t)) ~ Hb1 * _Ca_f(x,y,t) * (Bh - _CaB_hf(x,y,t)) - Hb2 * _CaB_hf(x,y,t),
	#These are non-linear voltage properties of the kv and ca channels
	Dt(mKV(x,y,t)) ~ αmKV(V(x,y,t)) * (1 - mKV(x,y,t)) - βmKV(V(x,y,t)) * mKV(x,y,t),
	Dt(hKV(x,y,t)) ~ αhKV(V(x,y,t)) * (1 - hKV(x,y,t)) - βhKV(V(x,y,t)) * hKV(x,y,t),
	Dt(mCa(x,y,t)) ~ αmCa(V(x,y,t)) * (1 - mCa(x,y,t)) - βmCa(V(x,y,t)) * mCa(x,y,t),
	#This is the non-linear term for the Calcium activate potassium
	Dt(mKCa(x,y,t)) ~ αmKCa(V(x,y,t)) * (1 - mKCa(x,y,t)) - βmKCa(V(x,y,t)) * mKCa(x,y,t),
	#Hyperpolarizing activated current
	Dt(C1(x,y,t)) ~ fhT(V(x,y,t), C1(x,y,t), C2(x,y,t), O1(x,y,t), O2(x,y,t), O3(x,y,t))[1], #This transforms the hMAT variable
	Dt(C2(x,y,t)) ~ fhT(V(x,y,t), C1(x,y,t), C2(x,y,t), O1(x,y,t), O2(x,y,t), O3(x,y,t))[2], #This transforms the hMAT variable
	Dt(O1(x,y,t)) ~ fhT(V(x,y,t), C1(x,y,t), C2(x,y,t), O1(x,y,t), O2(x,y,t), O3(x,y,t))[3], #This transforms the hMAT variable
	Dt(O2(x,y,t)) ~ fhT(V(x,y,t), C1(x,y,t), C2(x,y,t), O1(x,y,t), O2(x,y,t), O3(x,y,t))[4], #This transforms the hMAT variable
	Dt(O3(x,y,t)) ~ fhT(V(x,y,t), C1(x,y,t), C2(x,y,t), O1(x,y,t), O2(x,y,t), O3(x,y,t))[5], #This transforms the hMAT variable
	#Glutamate release
	#Dt(GL) ~ (-k_Ca * I_Ca(V) - GL) / τGL
]
# Set the boundary conditions and parameters ════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════#
bcs = [
	hv(x,y,tmin) ~  0.0      #Photons/μm
	V(x,y,tmin) ~  -35.0     #Membrane Voltage
	Dx(V(xmin, y, t)) ~ 0.0
	Dx(V(xmax, y, t)) ~ 0.0
	Dy(V(x, ymin, t)) ~ 0.0 #These boundary conditions matter the mose
	Dy(V(x, ymax, t)) ~ 0.0
	mKV(x,y,tmin) ~  0.430   #Rectifying potassium variable 1
	hKV(x,y,tmin) ~  0.999   #Rectifying potassium variable 2
	mCa(x,y,tmin) ~  0.436   #Voltage gated calcium current
	mKCa(x,y,tmin) ~  0.642  #Calcium activated potassium channel
	
	#Hyperpolarising currents
	C1(x,y,tmin) ~  0.646
	C2(x,y,tmin) ~  0.298
	O1(x,y,tmin) ~  0.0517
	O2(x,y,tmin) ~  0.00398
	O3(x,y,tmin) ~  0.000115

	#Phototransduction cascade
	R(x,y,tmin) ~  10000.0
	Ra(x,y,tmin) ~  0.0
	Trαβγ(x,y,tmin) ~  1000.0
	Trα(x,y,tmin) ~  0.0
	Trαi(x,y,tmin) ~  0.0
	Trβγ(x,y,tmin) ~  0.0
	PDE(x,y,tmin) ~  100.0
	PDEa(x,y,tmin) ~  0.0
	cGMP(x,y,tmin) ~  2.0
	GMP(x,y,tmin) ~  0.0

	Ca(x,y,tmin) ~  0.3 #Initial calcium
	CaB(x,y,tmin) ~  272.72727272727286 #Calcium Buffers

	#Calcium binding proteins
	_Ca_s(x,y,tmin) ~  0.0966 #Free calcium
	_Ca_f(x,y,tmin) ~  0.0966
	_CaB_ls(x,y,tmin) ~  80.929
	_CaB_hs(x,y,tmin) ~  29.068
	_CaB_lf(x,y,tmin) ~  80.929
	_CaB_hf(x,y,tmin) ~  29.068
	#GL(x,y,tmin) ~  0.0
];
domains = [
	x ∈ Interval(xmin, xmax)
	y ∈ Interval(ymin, ymax)
	t ∈ Interval(tmin, tmax)
]
dimensions = [x, y, t]
states = [
	hv(x,y,t),       #Photons/μm
	V(x,y,t),     #Membrane Voltage
	mKV(x,y,t),    #Rectifying potassium variable 1
	hKV(x,y,t),    #Rectifying potassium variable 2
	mCa(x,y,t),    #Voltage gated calcium current
	mKCa(x,y,t),   #Calcium activated potassium channel
	
	#Hyperpolarising currents
	C1(x,y,t), 
	C2(x,y,t), 
	O1(x,y,t), 
	O2(x,y,t), 
	O3(x,y,t), 

	#Phototransduction cascade
	R(x,y,t), 
	Ra(x,y,t), 
	Trαβγ(x,y,t), 
	Trα(x,y,t), 
	Trαi(x,y,t), 
	Trβγ(x,y,t), 
	PDE(x,y,t), 
	PDEa(x,y,t), 
	cGMP(x,y,t), 
	GMP(x,y,t), 

	Ca(x,y,t),  #Initial calcium
	CaB(x,y,t),  #Calcium Buffers

	#Calcium binding proteins
	_Ca_s(x,y,t),  #Free calcium
	_Ca_f(x,y,t), 
	_CaB_ls(x,y,t), 
	_CaB_hs(x,y,t), 
	_CaB_lf(x,y,t), 
	_CaB_hf(x,y,t), 
	#GL(x,y,tmin) =>  0.0
];
p0 = [
	αC => 0.6016,      #Collecting area    DEFAULT: 0.6016
	K1r => 50.0, #Regeneration of Rhodopsin
	K2f => 0.5, #Transducin activation by Rhodopsin
	K2r => 2.5, #Regeneration of Transducin
	K3f => 0.2,
	K3r => 5.0,
	K4 => 1.4,
	K5 => 0.8,
	gDARK => 5040.0,
	cm => 0.02,
	gL => 0.35,   #Leaky Conductance       DEFAULT:  0.3335 nS
	EL => -77.0, 
	gKV => 2.0,   #KV conductance          DEFAULT:  2.0 nS 
	EK => -74.0,  #Potassium Reversel      DEFAULT: -74.0 mV
	AMAX => 65.6,
	b => 0.25,
	γCa => 50.0,
	C0 => 0.1,      #Prevent calcium from falling below this limit
	er => 500,      #This is the limits of calcium buffering DEFAULT 500uM
	K6r => 0.2,     #Buffering of calcium    DEFAULT: 0.2 s^1 μM^-1 
	K6f => 0.8,     #Unbuffering of calcium  DEFAULT: 0.8 s^-1
	Kc => 0.1,      #Dimensionless reverse   DEFAULT: 0.1
	DCa => 6e-8,    #                        DEFAULT: 6.0e-8 dm^2 
	Lb1 => 0.4,     #                        DEFAULT: 0.4 s^-1
	Lb2 => 0.2,     #                        DEFAULT: 0.2 s^-1
	Hb1 => 100,     #                        DEFAULT: 100 s^-1 μM^-1
	Hb2 => 90,      #                        DEFAULT: 90 s^-1
	Bl => 500,      #                        DEFAULT: 500 μM
	Bh => 300,      #                        DEFAULT: 300 μm
	gCa => 0.7,     #Calicum maximal cond    DEFAULT: 0.7 nS
	ECa => -12.5,
	Cae => 0.01,    #                        DEFAULT: 0.01 μM
	J_ex => 20.0,   #                        DEFAULT: 20.0 pA
	J_ex2 => 20.0,  #                        DEFAULT: 20.0 pA
	K_ex => 2.3,    #                        DEFAULT: 2.3 μM 
	K_ex2 => 0.5,   #                        DEFAULT: 0.5 μM
	S1 => 3.142e-8, #const                   DEFAULT: 3.142e-8 dm^2
	δ => 3e-5,      #Calcium area            DEFAULT: 3.0e-5
	F => 9.648e4,   #This is Faraday const   DEFAULT: 9.648e4 C mol^-1
	V1 => 3.812e-13, #A second const		DEFAULT: 3.812e-13 dm^3
	V2 => 5.236e-13, #A third const		DEFAULT: 5.236e-13 dm^3
	gCl => 2.0,   #Chloride channel cond.  DEFAULT: 2.0 nS
	ECl => -20,   #Cl Reversal             DEFAULT: -20.0 mV 
	gKCa => 5.0,  #Calcium activated k     DEFAULT: 5.0 nS
	gH => 3.0,     #h channels conductance  DEFAULT: 3.0 nS
	EH => -32.0,   #h channel reversal      DEFAULT: -32.0 mV
	gGAP => 10.0
];
PC_PDE = PDESystem(PDEeqs, bcs, domains, dimensions, states, p0, name=:PC_SPDE)

# Discretize the model
print("Discretizing the model in: ")
discretization = MOLFiniteDifference([x => dx, y => dy], t)
GRID = get_discrete(PC_PDE, discretization) #Make a representation of the discrete map
@time prob = discretize(PC_PDE, discretization)

# Make the light stimulus ════════════════════════════════════════════════════════════════════════════════════════════════════════════════════════#
begin
	stim_start = 0.2
	stim_end = 0.4
	stimulus = zeros(11*11, length(tstops))
	stimulus[5*5, stim_start .< tstops .< stim_end] .= 0.0
	function affect!(integrator) 
		println("Simulating $(integrator.t)")
		integrator.u[1:11^2] .= stimulus[:, round(Int64, (integrator.t / dt) + 1)]
	end
	cbs = PresetTimeCallback(tstops, affect!)
end

# After loading the light stimulus, run the model
print("Running the model in: ")
@time sol = solve(prob, ROS3P(), callback = cbs, tstops = tstops)

#%% plot the solutions
fig, axs = plt.subplots(2)
TIME = sol.t
hvs = reshape(sol[hv(x,y,t)], 121, length(TIME))
vts = reshape(sol[V(x,y,t)], 121, length(TIME))

axs[1].plot(TIME, hvs')
axs[2].plot(TIME, vts')
axs[2].set_ylim()

# THis is my shelved PDE for the grid. It may be easier just to construct this as a ODE
#Load the dimensions at the top of the stack
dt = 0.01; dx = 0.1; dy = 0.1
xmin = ymin = tmin = 0.0
tmax = 100.0
xmax = ymax = 1.0
tspan = (tmin, tmax)

include("..\\src\\StarburstModel\\Variables.jl")
include("..\\src\\StarburstModel\\Domains.jl")
include("..\\src\\StarburstModel\\Equations.jl")

#set up the domains
@named SAC_PDE = PDESystem(SAC_eqs_PDE, SAC_Boundary, domains, dimensions, SAC_states_PDE, p0)
x_points = collect(xmin:dx:xmax)  # x grid points
y_points = collect(ymin:dy:ymax)  # y grid points

#randomize the points on the grid
#x_points = x_points .+ (rand(11)*0.01)
#x_points[1] = 0.0
#x_points[end] = 1.0
#y_points = y_points .+ (rand(11)*0.01)
#y_points[1] = 0.0
#y_points[end] = 1.0

discretization = MOLFiniteDifference([x => x_points, y => y_points], t)
GRID = get_discrete(SAC_PDE, discretization) #Make a representation of the discrete map
GRID[x]
@time SAC_prob_ODE = discretize(SAC_PDE, discretization);

discretization_EVEN = MOLFiniteDifference([x => dx, y => dy], t)
GRID = get_discrete(SAC_PDE, discretization_EVEN) #Make a representation of the discrete map
GRID[x]
@time SAC_prob_ODE = discretize(SAC_PDE, discretization_EVEN);

#=
sol = solve(SAC_prob_ODE)

nx = length(GRID[x])
ny = length(GRID[y])

SAC_prob_SPDE = SDEProblem(SAC_prob_ODE.f, SAC_noise_eqs_PDE, SAC_prob_ODE.u0, SAC_prob_ODE.tspan, SAC_prob_ODE.p)
sol = solve(SAC_prob_SPDE, SOSRI())
=#

using Revise
using ElectroPhysiology
using PhysiologyModeling

#3 Goals
#GOAL) Analysis of Physiological data
#GOAL) Speeding up the model
#GOAL) Analyzing the model
#GOAL) Testing models against physiology data

#%% GOAL: Testing models against physiology data
#I think a good next step is to add in the ability to open and compare data versus physiologyical data


#=[Solving a single SDE for tspan]==========================================#
#Extract and modify initial conditions
u0_dict = SAC_u0_dict()
u0 = extract_u0(u0_dict)

#Specify the timespan
tspan = (0.0, 300e3)

#Set the stimulus parameters
save_fn = "D:/Data/Analysis/2024_02_12_WT_Cell1_IC.png"
p0_dict = SAC_p0_dict() #Extract parameters
p0_dict["g_GABA"] = 0.0 #Null GABA conductance
p0_dict["g_ACh"] = 0.0 #Null ACh conductance
p0 = extract_p0(p0_dict)

#Set up the problem
prob = SDEProblem(SAC_ODE, noise1D, u0, tspan, p0)
@time sol = solve(prob, SOSRI(), reltol = 2e-2, abstol = 2e-2, progress = true, progress_steps = 1)

exp = Experiment(sol, dt = 1.0)
exp.t


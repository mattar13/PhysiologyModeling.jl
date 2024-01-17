#Activate the testing environment
using Pkg; Pkg.activate("test")

#I think a good next step is to add in the ability to open and compare data versus physiologyical data
# create a model output here
using Revise
using ElectroPhysiology, PhysiologyModeling

#First try opening some data
#Also we need to add some of the documentation to ElectroPhysiology.jl

#1) Run the model using the default settings
tspan = (0.0, 10e3)
p0 = extract_p0(SAC_p0_dict)
prob = SDEProblem(SAC_ODE, noise1D, vals_u0, tspan, p0)
@time sol = solve(prob, SOSRI(), reltol = 0.01, abstol = 0.01, progress = true, progress_steps = 1)

#Turn this into an experiment
sim_exp = Experiment(sol)


#Now we can open a experiment of our own
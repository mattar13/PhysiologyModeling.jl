using Revise
using ElectroPhysiology, PhysiologyModeling
using Pkg; Pkg.activate("test") #Activate the testing environment
using PhysiologyPlotting, GLMakie
using BenchmarkTools
using NLsolve
using BifurcationKit

#=============================================================================#
#%% Dimensional Analysis 
#    This is a 3D graph for which the phase planes are 2D. 
#    This shows the shape of the diffeq model
#    This section will result in a video showing the parameter Iapp change 
#=============================================================================#

tmin = 0.0
tmax = 300.0
dt = 0.1
p0_dict = SAC_p0_dict(keyset = :dynamical_analysis)
p0_dict["I_app"] = 0.0
p0 = extract_p0(p0_dict)
u0_dict = SAC_u0_dict()
u0 = extract_u0(u0_dict)
prob = ODEProblem(SAC_ODE, u0, (tmin, tmax), p0) #Create the problem

#The goal is to find the real numbers where fODE(u, p) = 0.0
function fODE(u, p) 
     du = similar(u)
     prob.f(du, u, p, 0.0)
     du
end

#Create the bifurcation problem from the fODE using the lens of the first parameter I_app
BifProb = BifurcationProblem(fODE, u0, p0, (@lens _[1]), 
     record_from_solution = (x, p) -> (v = x[2], n = x[3], iapp = p[1])
)

#Set the optiopns
p_min = -65.0
p_max = 50.0
opts_br = ContinuationPar(
     p_min = p_min, p_max = p_max, #This is the parameter range we want to observe
     max_steps = 10000, #Increase the number of steps does better at finding the unstable eqs
     detect_bifurcation = 3 #Detect bi
)
#Continuation acts like solve and solves the eq fODE(u, p) = 0.0
br = continuation(BifProb, PALC(), opts_br;	bothside = true, verbosity = 2)

#Extract the v, n, and i from the 
i_eq = map(i -> br.branch[i].iapp, 1:length(br))
v_eq = map(i -> br.branch[i].v, 1:length(br))
n_eq = map(i -> br.branch[i].n, 1:length(br))

#Run 200 simulations for the chosen parameters
n_traces = length(i_eq)
prng = i_eq
efunc(prob, i, repeat) = OneVarEnsembleProb(prob, i, repeat, prng)
ensemble_prob = EnsembleProblem(prob, prob_func = efunc)
@time sim = solve(ensemble_prob, EnsembleDistributed(), trajectories = n_traces, progress = true, progress_steps = 1, maxiters = 1e7);

#%%=============================================================================#
# Phase plane analysis  
#    Running a analysis of a single phase plane is a good place to start when running a dynamical analysis script    
#    This section will result in a video showing the parameter Iapp change 
#=============================================================================#

#Set the phase plane and nullcline parameters
nx = ny = 50 #Set the space over which to check the model
vmap = LinRange(-110.0, 10.0, nx) #Set the voltage map
nmap = LinRange(-0.1, 1.1, ny) #Set the nmap

#Contained within the ensemble simulation is the necessary models for phase plane analysis
sol = sim[1]
phase_map = phase_plane(sol.prob, vmap, nmap) #return the phase map
vt_vector = phase_map[:,:,1] #Extract the voltage vector
nt_vector = phase_map[:,:,2] #Extract the repol vector
strength = sqrt.(vt_vector .^ 2 .+ nt_vector .^ 2) #extract the strength
#Find the nullclines
vt_nc, nt_nc = find_nullclines(prob, vmap, nmap)
Time = sol.t[1]:dt:sol.t[end]
v_trace = map(t -> sol(t)[2], Time)
n_trace = map(t -> sol(t)[3], Time)

fp = br.branch[1]
v_br = fp.v
n_br = fp.n

#%%====================================Plot the solution====================================#
fig = Figure(size = (800, 800)); #Make the figure
g1 = fig[1,1] = GridLayout()
display(fig)

axA1 = Axis3(g1[1:2, 1], xlabel = "Iapp (pA)", ylabel = "V_eq (mV)", zlabel = "N_eq"); #Make axis 3
axA2 = Axis(g1[1,2], limits = (-65.0, 50.0, -110.0, 10.0), ylabel = "V_eq (mV)"); #Make axis 1
axA3 = Axis(g1[2,2], limits = (-65.0, 50.0, -0.1, 1.1), xlabel = "I applied (pA)", ylabel = "N_eq"); #Make axis 2

g2 = fig[2,1:2] = GridLayout()
axB1 = Axis(g2[1, 1], limits = (-110.0, 10.0, -0.1, 1.1), title = "Vt vector", xlabel = "Voltage (mV)", ylabel = "Repol (v_trace)")
axB2 = Axis(g2[1, 2], limits = (-110.0, 10.0, -0.1, 1.1), title = "Nt vector", xlabel = "Voltage (mV)", ylabel = "Repol (n_trace)")
axB3 = Axis(g2[1, 3], limits = (-110.0, 10.0, -0.1, 1.1), title = "Phase Plot", xlabel = "Voltage (mV)", ylabel = "Repol (n_trace)")

g3 = fig[3,1] = GridLayout()
axC1 = Axis(g3[1, 1], title = "I_clamp = 0.0pA", xlabel = "Time (ms)", ylabel = "Voltage (v_trace)")
axC2 = Axis(g3[2, 1], limits = (0.0, 300.0, -0.1, 1.1), xlabel = "Time (ms)", ylabel = "Repol (n_trace)")

g3b = fig[3,2] = GridLayout()
axC3 = Axis(g3b[1, 1], limits = (-110.0, 10.0, -0.1, 1.1), xlabel = "Voltage (mV)", ylabel = "Nt")

#Plot the bifurcation branches
lines!(axA1, i_eq, v_eq, n_eq, linewidth = 1.0, color = :black)
lines!(axA2, i_eq, v_eq, linewdith = 1.0, color = :black)
lines!(axA3, i_eq, n_eq, linewidth = 1.0, color = :black)

#Plot the branch points
for i in br.specialpoint
     label = "I_$(String(i.type)) = $(round(i.param, digits = 2)) pA"
     scatter!(axA1, i.param, i.x[2], i.x[3])
     scatter!(axA2, i.param, i.x[2], label = label)
     scatter!(axA3, i.param, i.x[3])
end
fig[1,2] = Legend(fig, ax2, "Branchpoints", framevisible = false)

vfield1 = heatmap!(axB1, vmap, nmap, vt_vector, alpha = 0.5)
nfield2 = heatmap!(axB2, vmap, nmap, nt_vector, alpha = 0.5)
phase3 = heatmap!(axB3, vmap, nmap, strength, alpha = 0.5)

nNC1 = lines!(axB1, vmap, nt_nc, color = :red, linewidth = 2.0)
vNC2 = lines!(axB2, vmap, vt_nc, color = :blue, linewidth = 2.0)
nNC_PhaseNt = lines!(axB3, vmap, nt_nc, color = :red, linewidth = 2.0)
vNC_PhaseVt = lines!(axB3, vmap, vt_nc, color = :blue, linewidth = 2.0)
phase_trace = lines!(axB3, v_trace, n_trace, color = :black, linewidth = 2.0, linestyle = :dashdot)

fp1 = scatter!(axB3, v_eq[1], n_eq[1], markersize = 10.0)
lTR1 = lines!(axC1, Time, v_trace, color = :red, linewidth = 2.0)
la2 = lines!(axC2, Time, n_trace, color = :blue, linewidth = 2.0)
lb1 = lines!(axC3, v_trace, n_trace, color = :black, linewidth = 2.0, linestyle = :dashdot)


#%% Animate this section of the model
record(fig, "data/DynamicAnalysis.mp4", enumerate(sim), framerate = 10) do (i, sol)
     println(i)

     phase_map = phase_plane(sol.prob, vmap, nmap) #return the phase map
     vt_vector = phase_map[:,:,1] #Extract the voltage vector
     nt_vector = phase_map[:,:,2] #Extract the repol vector
     strength = sqrt.(vt_vector .^ 2 .+ nt_vector .^ 2) #extract the strength
     #Find the nullclines
     vt_nc, nt_nc = find_nullclines(sol.prob, vmap, nmap)
     Time = sol.t[1]:dt:sol.t[end]
     v_trace = map(t -> sol(t)[2], Time)
     n_trace = map(t -> sol(t)[3], Time)

     fp = br.branch[i]
     #= If there is a certain criteria, we want to plot this section
     lines!(axA1, [prng[i]], v_trace, n_trace, color = :red, alpha = 0.5, depth_shhift = 0)
     lines!(axA2, [prng[i]], v_trace, color = :red, alpha = 0.25, depth_shhift = 0)
     lines!(axA3, [prng[i]], n_trace, color = :red, alpha = 0.25, depth_shhift = 0)
     =#
     axC1.title = "I_clamp = $(round(prng[i], digits = 2)) pA"

     vfield1[3] = vt_vector
     nfield2[3] = nt_vector
     phase3[3] = strength

     nNC1[2] = nt_nc;
     vNC2[2] = vt_nc;
     nNC_PhaseNt[2] = nt_nc
     vNC_PhaseVt[2] = vt_nc
     phase_trace[1] = v_trace; phase_trace[2] = n_trace


end


#=============================================================================#
#%% Codimensional Analysis 
#    This is a 3D graph for which the phase planes are 2D. 
#    This shows the shape of the diffeq model
#    This section will result in a video showing the parameter Iapp change 
#=============================================================================#

tmin = 0.0
tmax = 300.0
dt = 0.1
p0_dict = SAC_p0_dict(keyset = :dynamical_analysis)
p0_dict["I_app"] = 0.0
p0 = extract_p0(p0_dict)
u0_dict = SAC_u0_dict()
u0 = extract_u0(u0_dict)
prob = ODEProblem(SAC_ODE, u0, (tmin, tmax), p0) #Create the problem

#The goal is to find the real numbers where fODE(u, p) = 0.0
function fODE(u, p) 
     du = similar(u)
     prob.f(du, u, p, 0.0)
     du
end

#Set the optiopns
p_min = -65.0
p_max = 50.0
opts_br = ContinuationPar(
     p_min = p_min, p_max = p_max, #This is the parameter range we want to observe
     max_steps = 10000, #Increase the number of steps does better at finding the unstable eqs
     detect_bifurcation = 3 #Detect bi
)
#Continuation acts like solve and solves the eq fODE(u, p) = 0.0
br = continuation(BifProb, PALC(), opts_br;	bothside = true, verbosity = 2)

#Extract the v, n, and i from the 
i_eq = map(i -> br.branch[i].iapp, 1:length(br))
v_eq = map(i -> br.branch[i].v, 1:length(br))

#====================================Plot the solution====================================#
fig = Figure(size = (800, 400)); #Make the figure
ax1 = Axis(fig[1,1], ylabel = "v_trace (mV)"); #Make axis 1
#Plot the bifurcation branches
lines!(ax1, i_eq, v_eq, linewdith = 1.0, color = :black)

#Plot the 
for i in br.specialpoint
     label = "I_$(String(i.type)) = $(round(i.param, digits = 2)) pA"
     scatter!(ax1, i.param, i.x[2], label = label)
end
fig[1,2] = Legend(fig, ax2, "Branchpoints", framevisible = false)

Time = 200.0:dt:tmax
for (i, sol) in enumerate(sim)
     v_trace = map(t -> sol(t)[2], Time)
     n_trace = map(t -> sol(t)[3], Time)
     lines!(ax1, [prng[i]], v_trace, color = :red, alpha = 0.25, depth_shhift = 0)
end
display(fig)
using Revise
using Pkg; Pkg.activate(".")
using ElectroPhysiology, PhysiologyModeling
Pkg.activate("test") #Activate the testing environment
using PhysiologyPlotting, GLMakie 
using BifurcationKit

#=============================================================================#
#%% Setting up the experiment with a dynamical system
#    This is a 3D graph for which the phase planes are 2D. 
#    This shows the shape of the diffeq model
#    This section will result in a video showing the parameter Iapp change 
#=============================================================================#
tmin = 0.0 #Set the minimum time for the trajectory
tmax = 300.0 #Set the maximum time for the trajectory
dt = 0.1
p0_dict = SAC_p0_dict(keyset = :dynamical_analysis)
p0_dict["I_app"] = 0.0
p0 = extract_p0(p0_dict)
u0_dict = SAC_u0_dict()
u0 = extract_u0(u0_dict, mode = :DynamicalAnalysis)
prob = ODEProblem(DynamicSAC, u0, (tmin, tmax), p0) #Create the problem

#The goal is to find the real numbers where fODE(u, p) = 0.0
function fODE(u, p) 
     du = similar(u)
     prob.f(du, u, p, 0.0)
     du
end

#=============================================================================#
# Dimensional Analysis 
#    This plots all the equilibria over a changing model
#    This shows the shape of the diffeq model
#    This section will result in a video showing the parameter Iapp change 
#=============================================================================#
#Create the bifurcation problem from the fODE using the lens of the first parameter I_app
BifProb = BifurcationProblem(fODE, u0, p0, (@lens _[1]), 
     record_from_solution = (x, p) -> (v = x[1], n = x[2], iapp = p[1])
)

#Set the optiopns
p_min = -65.0
p_max = 200.0
opts_br = ContinuationPar(
     p_min = p_min, p_max = p_max, #This is the parameter range we want to observe
     max_steps = 10000, #Increase the number of steps does better at finding the unstable eqs
     detect_bifurcation = 3 #Detect bi
)
#Continuation acts like solve and solves the eq fODE(u, p) = 0.0
br = continuation(BifProb, PALC(), opts_br; bothside = true, verbosity = 2)

#Extract the v, n, and i from the continuation
i_eq = map(i -> br.branch[i].iapp, 1:length(br))
v_eq = map(i -> br.branch[i].v, 1:length(br))
n_eq = map(i -> br.branch[i].n, 1:length(br))

#=============================================================================#
# Multiple simulation range 
#    Running a analysis of a single phase plane is a good place to start when running a dynamical analysis script    
#    This section will result in a video showing the parameter Iapp change 
#=============================================================================#

#Run 200 simulations for the chosen parameters
n_traces = 200
prng = LinRange(p_min, p_max, n_traces)
efunc(prob, i, repeat) = OneVarEnsembleProb(prob, i, repeat, prng)
ensemble_prob = EnsembleProblem(prob, prob_func = efunc)
@time sim = solve(ensemble_prob, EnsembleDistributed(), trajectories = n_traces, progress = true, progress_steps = 1, maxiters = 1e7);

#=============================================================================#
# Phase plane analysis over the simulation range
#    Running a analysis of a single phase plane is a good place to start when running a dynamical analysis script    
#    This section will result in a video showing the parameter Iapp change 
#=============================================================================#
#Set the phase plane and nullcline parameters
nx = ny = 50 #Set the space over which to check the model
vmap = LinRange(-110.0, 10.0, nx) #Set the voltage map
nmap = LinRange(-0.1, 1.1, ny) #Set the nmap

#for each solution in the simulation calculate a phase plane, nullcline, and fixedpoints
v_traces = []; n_traces = []
vt_vectors = []; nt_vectors = []
v_fps = []; n_fps = []
v_ncs = []; n_ncs = []
Time = tmin:dt:tmax
for (idx, sol) in enumerate(sim)
     
     v_trace = map(t-> sol(t)[1], Time)
     n_trace = map(t-> sol(t)[2], Time)
     push!(v_traces, v_trace)
     push!(n_traces, n_trace)

     prob_i = sol.prob
     phase_map = phase_plane(prob_i, vmap, nmap) #return the phase map
     push!(vt_vectors, phase_map[:,:,1]) #Extract the voltage vector
     push!(nt_vectors, phase_map[:,:,2]) #Extract the repol vector

     ncs = find_nullclines(prob_i, vmap, nmap) #find the nullclines
     push!(v_ncs, ncs[1])
     push!(n_ncs, ncs[2])

     fps = find_fixed_points(prob_i, vmap, nmap) #find the fixed points
     push!(v_fps, map(FP -> FP[1], fps))
end
strengths = map(i -> sqrt.(vt_vectors[i] .^ 2 .+ nt_vectors[i] .^ 2), eachindex(vt_vectors))
v_fps = map(i -> !isempty(i) ? i : [NaN], v_fps)
n_fps = map(i -> !isempty(i) ? i : [NaN], n_fps)

#%%====================================Plot the solution====================================#
fig = Figure(size = (800, 800)); #Make the figure
g1 = fig[1,1] = GridLayout()
display(fig)

#Make the layout of the plot
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

#Plot the bifurcation branches first
lines!(axA1, i_eq, v_eq, n_eq, linewidth = 1.0, color = :black)
lines!(axA2, i_eq, v_eq, linewdith = 1.0, color = :black)
lines!(axA3, i_eq, n_eq, linewidth = 1.0, color = :black)

#Plot the branching points
for i in br.specialpoint
     label = "I_$(String(i.type)) = $(round(i.param, digits = 2)) pA"
     scatter!(axA1, i.param, i.x[1], i.x[2])
     scatter!(axA2, i.param, i.x[1], label = label)
     scatter!(axA3, i.param, i.x[2])
end
fig[1,2] = Legend(fig, axA2, "Branchpoints", framevisible = false)

#Plot the vector fields of Vt and Nt
vfield1 = heatmap!(axB1, vmap, nmap, vt_vectors[1], alpha = 0.5)
nfield2 = heatmap!(axB2, vmap, nmap, nt_vectors[1], alpha = 0.5)
phase3 = heatmap!(axB3, vmap, nmap, strengths[1], alpha = 0.5)

#Plot the nullclines
nNC1 = lines!(axB1, vmap, n_ncs[1], color = :red, linewidth = 2.0)
vNC2 = lines!(axB2, vmap, v_ncs[1], color = :blue, linewidth = 2.0)
nNC_PhaseNt = lines!(axB3, vmap, n_ncs[1], color = :red, linewidth = 2.0)
vNC_PhaseVt = lines!(axB3, vmap, v_ncs[1], color = :blue, linewidth = 2.0)
phase_trace = lines!(axB3, v_traces[1], n_traces[1], color = :black, linewidth = 2.0, linestyle = :dashdot)

#Plot the fixed points
points = Point2f.(v_fps[1], n_fps[1])
fp1 = scatter!(axB3, points, markersize = 10.0)
la1 = lines!(axC1, Time, v_traces[1], color = :red, linewidth = 2.0)
la2 = lines!(axC2, Time, n_traces[1], color = :blue, linewidth = 2.0)
lb1 = lines!(axC3, v_traces[1], n_traces[1], color = :black, linewidth = 2.0, linestyle = :dashdot)

#Animate this section of the model
x_rng = 1:length(sim)

record(fig, "test/SAC_model_tests/data/DynamicAnalysis.mp4", x_rng, framerate = 10) do i
     println(i)
     sol = sim[i]
     axC1.title = "I_clamp = $(round(i_eq[i], digits = 2)) pA"

     vfield1[3] = vt_vectors[i]
     nfield2[3] = nt_vectors[i]
     phase3[3] = strengths[i]

     nNC1[2] = n_ncs[i];
     vNC2[2] = v_ncs[i];
     nNC_PhaseNt[2] = n_ncs[i]
     vNC_PhaseVt[2] = v_ncs[i]
     phase_trace[1] = v_traces[i]; 
     phase_trace[2] = n_traces[i]
     
     points = Point2f.(v_fps[i], n_fps[i])
     fp1[1] = points
     la1[2] = v_traces[i]
     la2[2] = n_traces[i]
     lb1[1] = v_traces[i]; 
     lb1[2] = n_traces[i]
end

#=============================================================================#
#%% Codimensional Analysis 
#    This is a 3D graph for which the phase planes are 2D. 
#    This shows the shape of the diffeq model
#    This section will result in a video showing the parameter Iapp change 
#=============================================================================#

#Create the bifurcation problem from the fODE using the lens of the first parameter I_app
BifProb = BifurcationProblem(fODE, u0, p0, (@lens _[1]), 
     record_from_solution = (x, p) -> (v = x[1], n = x[2], iapp = p[1])
)
 
#Set the options
p_min = -65.0
p_max = 200.0
opts_br = ContinuationPar(
     p_min = p_min, p_max = p_max, #This is the parameter range we want to observe
     max_steps = 10000, #Increase the number of steps does better at finding the unstable eqs
     detect_bifurcation = 3 #Detect bi
)
#Continuation acts like solve and solves the eq fODE(u, p) = 0.0
br = continuation(BifProb, PALC(), opts_br; bothside = true)
br
#Extract the v, n, and i from the continuation
i_eq = map(i -> br.branch[i].iapp, 1:length(br))
v_eq = map(i -> br.branch[i].v, 1:length(br))
n_eq = map(i -> br.branch[i].n, 1:length(br))

#%% Do the codim 2 analysis
p2_min = 5.0
p2_max = 200.0
c2_opts = ContinuationPar(
     p_min = p2_min, p_max = p2_max, #This is the parameter range we want to observe
     detect_bifurcation = 2, #Detect bi
     n_inversion = 10, max_steps = 10000, max_bisection_steps = 55
)
cont_par = par_idx("g_Ca")
sn_codim2 = continuation(br, 2, (@lens _[cont_par]), c2_opts,
     update_minaug_every_step = 1, 
     detect_codim2_bifurcation = 2, 
     record_from_solution = (x, p) -> (v = x[1], n = x[2], gc = p[1])
)
#Extract the v, n, and i from the continuation
sn_codim2.branch
i_eq_c2 = map(i -> sn_codim2.branch[i].p1, 1:length(sn_codim2))
gc_eq_c2 = map(i -> sn_codim2.branch[i].p2, 1:length(sn_codim2))
v_eq_c2 = map(i -> sn_codim2.branch[i].v, 1:length(sn_codim2))
n_eq_c2 = map(i -> sn_codim2.branch[i].n, 1:length(sn_codim2))

#%%====================================Plot the solution====================================#
fig = Figure(size = (750, 750)); #Make the figure
#g1 = fig[1,1] = GridLayout()
display(fig)

#Make the layout of the plot
axA1 = Axis3(fig[1, 1], xlabel = "Iapp (pA)", ylabel = "V_eq (mV)", zlabel = "N_eq"); #Make axis 3
axA2 = Axis(fig[2,1:2], #=limits = (-65.0, 50.0, -110.0, 10.0),=# ylabel = "V_eq (mV)"); #Make axis 1
axA3 = Axis(fig[3,1:2], #=limits = (-65.0, 50.0, -0.1, 1.1),=# xlabel = "I applied (pA)", ylabel = "N_eq"); #Make axis 2
axA4 = Axis(fig[4,1:2], #=limits = (-65.0, 50.0, -0.1, 1.1),=# xlabel = "I applied (pA)", ylabel = "N_eq"); #Make axis 2

#Plot the bifurcation branches first
lines!(axA1, i_eq, v_eq, n_eq, linewidth = 1.0, color = :black)
lines!(axA2, i_eq, v_eq, linewdith = 1.0, color = :black)
lines!(axA3, i_eq, n_eq, linewidth = 1.0, color = :black)

lines!(axA1, i_eq_c2, v_eq_c2, n_eq_c2, linewidth = 1.0, color = :red)
lines!(axA2, i_eq_c2, v_eq_c2, linewdith = 1.0, color = :red)
lines!(axA3, i_eq_c2, n_eq_c2, linewidth = 1.0, color = :red)
lines!(axA4, gc_eq_c2, i_eq_c2, linewidth = 1.0, color = :red)

#Plot the branching points
for i in br.specialpoint
     label = "I_$(String(i.type)) = $(round(i.param, digits = 2)) pA"
     scatter!(axA1, i.param, i.x[1], i.x[2])
     scatter!(axA2, i.param, i.x[1], label = label)
     scatter!(axA3, i.param, i.x[2])
end
fig[1,2] = Legend(fig, axA2, "Branchpoints", framevisible = false)

display(fig)
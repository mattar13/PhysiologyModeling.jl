using Revise
using ElectroPhysiology, PhysiologyModeling
using Pkg; Pkg.activate("test") #Activate the testing environment
using PhysiologyPlotting, GLMakie
using BenchmarkTools
using NLsolve

#=============================================================================#
#%% Phase plane analysis  
#    Running a analysis of a single phase plane is a good place to start when running a dynamical analysis script    
#    This section will result in a video showing the parameter Iapp change 
#=============================================================================#

#Set up the figure first
fig2 = Figure(size = (1200, 800))
ga = fig2[1, 1:2] = GridLayout() #The first grid is the two voltage solutions
axA1 = Axis(ga[1, 1], limits = (0.0, 300.0, -90.0, 20.0), title = "I_clamp = 0.0pA", xlabel = "Time (ms)", ylabel = "Voltage (Vt)")
axA2 = Axis(ga[2, 1], limits = (0.0, 300.0, -0.1, 1.1), xlabel = "Time (ms)", ylabel = "Repol (Nt)")
gb = fig2[1,3] = GridLayout()
axB1 = Axis(gb[1, 1], limits = (-90.0, 20.0, -0.1, 1.1), xlabel = "Voltage (mV)", ylabel = "Repol (Nt)")
# Now plot the phase diagrams
gC = fig2[2, 1:3] = GridLayout()
axC1 = Axis(gC[1, 1], limits = (-90.0, 20.0, -0.1, 2.1), title = "Vt vector", xlabel = "Voltage (mV)", ylabel = "Repol (Vt)")
axC2 = Axis(gC[1, 2], limits = (-90.0, 20.0, -0.1, 2.1), title = "Nt vector", xlabel = "Voltage (mV)", ylabel = "Repol (Nt)")
axC3 = Axis(gC[1, 3], limits = (-90.0, 20.0, -0.1, 2.1), title = "Phase Plot", xlabel = "Voltage (mV)", ylabel = "Repol (Nt)")

#Load the model problem ics and parameters
u0_dict = SAC_u0_dict() #Load parameters
p0_dict = SAC_p0_dict(keyset = :dynamical_analysis) #Load parameters
p0_dict["I_app"] = 10.0 #Set the I_app
u0 = extract_u0(u0_dict) #extract the initial conditions
p0 = extract_p0(p0_dict) #Extract the initial parameters
prob = ODEProblem(SAC_ODE, u0, (0.0, 300.0), p0) #Create the problem
sol = solve(prob) #Solve the problem

#Extract the time and solution
dt = 0.001
Time = sol.t[1]:dt:sol.t[end]
vt = map(t -> sol(t)[2], Time)
nt = map(t -> sol(t)[3], Time)

#Set the phase plane and nullcline parameters
nx = ny = 50 #Set the space over which to check the model
vmap = LinRange(-90.0, 20.0, nx) #Set the voltage map
nmap = LinRange(-2.0, 2.0, ny) #Set the nmap
phase_map = phase_plane(prob, vmap, nmap) #return the phase map
vt_vector = phase_map[:,:,1] #Extract the voltage vector
nt_vector = phase_map[:,:,2] #Extract the repol vector
strength = sqrt.(vt_vector .^ 2 .+ nt_vector .^ 2) #extract the strength

#Find the nullclines and fixed points
vt_nc, nt_nc = find_nullclines(prob, vmap, nmap)
fixed_points = find_fixed_points(prob, vmap, nmap, precision = 5)
map(fp -> equilibrium_stability(prob, fp), fixed_points)
fixed_points_pt = Point2f.(map(pt -> (pt[2], pt[3]), fixed_points))

la1 = lines!(axA1, Time, vt, color = :red, linewidth = 2.0)
la2 = lines!(axA2, Time, nt, color = :blue, linewidth = 2.0)
lb1 = lines!(axB1, vt, nt, color = :black, linewidth = 2.0, linestyle = :dashdot)

hm1 = heatmap!(axC1, vmap, nmap, vt_vector, alpha = 0.5)
hm2 = heatmap!(axC2, vmap, nmap, nt_vector, alpha = 0.5)
hm3 = heatmap!(axC3, vmap, nmap, strength, alpha = 0.5)

lc1 = lines!(axC1, vmap, nt_nc, color = :red, linewidth = 2.0)
lc2 = lines!(axC2, vmap, vt_nc, color = :blue, linewidth = 2.0)
lcPhaseNt = lines!(axC3, vmap, nt_nc, color = :red, linewidth = 2.0)
lcPhaseVt =lines!(axC3, vmap, vt_nc, color = :blue, linewidth = 2.0)
lcTrace = lines!(axC3, vt, nt, color = :black, linewidth = 2.0, linestyle = :dashdot)

#Plot the fixed points
fp1 = scatter!(axC3, fixed_points_pt)

#Animate the solution
xrng = LinRange(-20.0, 100.0, 100)
record(fig2, "test/SAC_model_tests/PhasePlane.mp4", xrng, framerate = 10) do x
     println(x)
     p0_dict = SAC_p0_dict(keyset = :dynamical_analysis) #Load parameters
     p0_dict["I_app"] = x #Set the I_app
     u0 = extract_u0(u0_dict) #extract the initial conditions
     p0 = extract_p0(p0_dict) #Extract the initial parameters
     prob = ODEProblem(SAC_ODE, u0, (0.0, 300.0), p0) #Create the problem
     sol = solve(prob) #Solve the problem
     #Extract the time and solution
     dt = 0.001
     Time = sol.t[1]:dt:sol.t[end]
     vt = map(t -> sol(t)[2], Time)
     nt = map(t -> sol(t)[3], Time)
     #Set the phase plane and nullcline parameters
     phase_map = phase_plane(prob, vmap, nmap) #return the phase map
     vt_vector = phase_map[:,:,1] #Extract the voltage vector
     nt_vector = phase_map[:,:,2] #Extract the repol vector
     strength = sqrt.(vt_vector .^ 2 .+ nt_vector .^ 2) #extract the strength

     # finding nullclines
     vt_nc, nt_nc = find_nullclines(prob, vmap, nmap)
     fixed_points = find_fixed_points(prob, vmap, nmap)

     axA1.title = "I_clamp = $(round(x, digits = 2)) pA"
     la1[2] = vt
     la2[2] = nt
     lb1[1] = vt; lb1[2] = nt

     hm1[3] = vt_vector
     hm2[3] = nt_vector
     hm3[3] = strength

     lc1[2] = nt_nc;
     lc2[2] = vt_nc;
     lcPhaseNt[2] = nt_nc
     lcPhaseVt[2] = vt_nc
     lcTrace[1] = vt; lcTrace[2] = nt
     fp1[1] = Point2f.(fixed_points)
end

#=============================================================================#
#%% Dimensional Analysis 
#    This is a 3D graph for which the phase planes are 2D. 
#    This shows the shape of the diffeq model
#    This section will result in a video showing the parameter Iapp change 
#=============================================================================#
fig1 = Figure(size = (1600, 800))
ga = fig1[1, 1] = GridLayout()
axA1 = Axis(ga[1, 1], xlabel = "Time (ms)", ylabel = "Voltage (Vt)")
axA2 = Axis(ga[2, 1], xlabel = "Time (ms)", ylabel = "Repol (Nt)")

gb = fig1[1, 2] = GridLayout()
axB1 = Axis3(gb[1,1], title = "K Repolarization (Nt)", xlabel = "I_app", ylabel = "Voltage. (Vt) ", zlabel = "K Repol. (N)")

# Set some initial parameters
u0_dict = SAC_u0_dict()
p0_dict = SAC_p0_dict(keyset = :dynamical_analysis)
u0_dict["v"] = -65.0

# After setting the inital parameters, extract them
u0 = extract_u0(u0_dict)
p0 = extract_p0(p0_dict)
prob = ODEProblem(SAC_ODE, u0, (0.0, 300.0), p0) #Create the problem

#Set up the ensemble solution
n_traces = 100 #Number of traces
p_rng = LinRange(-20, 250, n_traces) #Set the parameter range
efunc(prob, i, repeat) = OneVarEnsembleProb(prob, i, repeat, p_rng)
ensemble_prob = EnsembleProblem(prob, prob_func = efunc)

#Solve the ensemble solution
@time sim = solve(ensemble_prob, EnsembleDistributed(), trajectories = n_traces, progress = true, progress_steps = 1, maxiters = 1e7);

# Plot the results
for (i, sol) in enumerate(sim)
     println(i)
     x = p_rng[i]
     dt = 0.01
     Time = sol.t[1]:dt:sol.t[end]
     vt = map(t -> sol(t)[2], Time)
     nt = map(t -> sol(t)[3], Time)
     ct = map(t -> sol(t)[6], Time)

     lines!(axA1, Time, vt, colormap = :rainbow, color = [x], colorrange = (p_rng[1], p_rng[end]))
     lines!(axA2, Time, nt, colormap = :rainbow, color = [x], colorrange = (p_rng[1], p_rng[end]))
     lines!(axB1, [x], vt, nt, colormap = :rainbow, color = [x], colorrange = (p_rng[1], p_rng[end]))
end
display(fig1)
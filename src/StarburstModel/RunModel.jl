function runODEModel(tmax; alg = Tsit5())
     #Load the dimensions at the top of the stack
     tspan = (0.0, tmax)
     @named SAC_ODE = ODESystem(SAC_eqs_ODE, t, SAC_states_ODE, SAC_parameters, tspan = tspan)
     prob = ODEProblem(SAC_ODE)
     sol = solve(prob, alg)
     #Load the solution into a Experiment
     Experiment(sol)
end


function runSDEModel(tmax)
     #Load the dimensions at the top of the stack
     tspan = (0.0, tmax)
     @named SAC_ODE = ODESystem(SAC_eqs_ODE, t, SAC_states_ODE, SAC_parameters)
     @named SAC_SDE = SDESystem(SAC_ODE, SAC_noise_eqs, tspan = tspan)
     prob = SDEProblem(SAC_SDE)
     sol = solve(prob, SOSRI())
     #Make an updated Experiment
     Experiment(sol)
end

function runPDEModel(tmax, xmax, ymax; dt = 0.01, dx = 0.1, dy = 0.1)
     xmin = ymin = tmin = 0.0
     tmax = 100.0
     xmax = ymax = 1.0
     tspan = (tmin, tmax)
     #Set the boundary conditions
     SAC_Boundary = [
          Î_ext(x, y, tmin) ~ 0.0
          v̂(x, y, tmin) ~ -63.6 
          n̂(x, y, tmin) ~ 0.000 
          m̂(x, y, tmin) ~ 0.062 
          ĥ(x, y, tmin) ~ 0.550 
          ĉ(x, y, tmin) ~ 0.085 
          â(x, y, tmin) ~ 0.026 
          b̂(x, y, tmin) ~ 0.000 
          ê(x, y, tmin) ~ 0.066 
          î(x, y, tmin) ~ 0.053     
          Ŵ(x, y, tmin) ~ 0.000 #Initial conditions
     
          Dx(ê(xmin, y, t)) ~ 0.0
          Dx(ê(xmax, y, t)) ~ 0.0
          Dy(ê(x, ymin, t)) ~ 0.0
          Dy(ê(x, ymax, t)) ~ 0.0

          Dx(î(xmin, y, t)) ~ 0.0
          Dx(î(xmax, y, t)) ~ 0.0
          Dy(î(x, ymin, t)) ~ 0.0     
          Dy(î(x, ymax, t)) ~ 0.0
     ]

     domains = [
          x ∈ IntervalDomain(xmin, xmax)
          y ∈ IntervalDomain(ymin, ymax)
          t ∈ IntervalDomain(tmin, tmax)
     ]
     @named SAC_PDE = PDESystem(SAC_eqs_PDE, SAC_Boundary, domains, dimensions, SAC_states_PDE, p0)
     discretization = MOLFiniteDifference([x => dx, y => dy], t)
     GRID = get_discrete(SAC_PDE, discretization) #Make a representation of the discrete map
     @time prob = discretize(SAC_PDE, discretization);
     sol = solve(prob)
     return sol     
end
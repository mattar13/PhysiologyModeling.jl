using Pkg; Pkg.activate(".")
using PhysiologyModeling
using FFTW

# ---------- load packages ----------
include("auxillary_functions.jl")
include("models.jl")
include("grid_discretization.jl")
include("parameters.jl")

# ---------- load plotting packages ----------
using Pkg; Pkg.activate("test")
using GLMakie
using Statistics
#%% -Test the grid-based model with a single release site in the center.

# Grid parameters
nx, ny = 50, 50
dx, dy = 1.0, 1.0
D = 0.1  # Diffusion coefficient
release_sites = [(25, 25)]  # Center of the grid
grid_params = GridParameters(nx, ny, dx, dy, D, release_sites)
grid_params.fdm_operator

#Create the parameters for the dopamine grid
p_grid = (p.krel, p.kclear)

#Create the initial conditions for the dopamine grid
u0 = zeros(nx*ny)
# Calculate middle point indices
mid_x = div(nx, 2) + 1  # Integer division by 2, then add 1
mid_y = div(ny, 2) + 1
# Convert 2D index to 1D index
mid_idx = (mid_y - 1) * nx + mid_x
u0[mid_idx] = 1.0

tspan = (0.0, 1000.0)

prob = ODEProblem(
    (du, u, p, t) -> update_dopamine_grid!(du, u, p, t, grid_params),
    u0, tspan, p_grid
)

sol = solve(prob, Tsit5(), saveat=1.0)  # Save more frequently for smoother animation

# Create 3D array of solutions
DA_grid = cat(map(t -> reshape(sol(t), nx, ny), sol.t)...; dims=3)

#%% Create animation
fig = Figure(size=(800, 800))
ax = Axis(fig[1, 1], title="Dopamine Concentration Evolution", aspect=1.0)

# Calculate the maximum value across all frames for consistent color limits
max_val = maximum(DA_grid)
hm = heatmap!(ax, DA_grid[:,:,1], colormap=:viridis, colorrange=(0, max_val))
Colorbar(fig[1, 2], hm)

# Add release site markers
for (i, j) in release_sites
    scatter!(ax, [i], [j], color=:red, markersize=15, label="Release Site")
end

# Set axis limits and labels
xlims!(ax, 0.5, nx+0.5)
ylims!(ax, 0.5, ny+0.5)
axislegend(ax)

# Create animation
framerate = 30
record(fig, "dopamine_evolution.mp4", 1:size(DA_grid, 3); framerate=framerate) do frame
    println(frame)
    hm[1] = DA_grid[:,:,frame]  # Update the heatmap data
    ax.title = "Dopamine Concentration at t = $(round(sol.t[frame], digits=1))"
end

#%% Run the model with the dopamine grid Initial conditions
n_grid = nx * ny
u0 = zeros(2 + n_grid + 4)  # V, Ca, DA_grid, P, Gi, cAMP, PKA
u0[1] = -65.0  # Initial voltage
u0[2] = 0.0    # Initial calcium

# Time span
tspan = (0.0, 1000.0)

# Create the ODE problem
prob = ODEProblem((du, u, p, t) -> dopaminergic_autoreceptor_grid!(du, u, p, t, grid_params, :fdm),
                    u0, tspan, p)

# Solve
sol = solve(prob, Tsit5(), saveat=10.0)

# Plot results using GLMakie
fig = Figure(size=(1200, 400))
times = [100, 500, 1000]  # Time points to plot

for (i, t) in enumerate(times)
    idx = findfirst(x -> x >= t, sol.t)
    DA_grid = reshape(sol.u[idx][3:end-4], nx, ny)
    
    ax = Axis(fig[1, i], title="t = $t", aspect = :equal)
    heatmap!(ax, DA_grid, colormap=:viridis)
    Colorbar(fig[1, i+1], limits=(0, maximum(DA_grid)))
end

display(fig)

#%% ---------- test multiple release sites ----------
"""
    test_multiple_release_sites()
Test the grid-based model with multiple release sites.
"""
function test_multiple_release_sites()
    # Grid parameters
    nx, ny = 50, 50
    dx, dy = 1.0, 1.0
    D = 0.1
    release_sites = [(15,15), (35,15), (25,35)]  # Three release sites
    
    grid_params = GridParameters(nx, ny, dx, dy, D, release_sites)
    
    # Use same parameters as single release site test
    p = (
        Cm=1.0, gL=0.1, EL=-65.0, gCa=1.0, ECa=120.0,
        Vhalf=-20.0, kV=10.0, τCa=100.0, αCa=0.1,
        krel=1.0, kclear=0.1, kon=1.0, koff=0.1,
        kG=1.0, kGdeg=0.1, G50=1.0, nGi=2.0,
        kACbasal=1.0, kcAMPdeg=0.1, kPKA=1.0, kPKAdeg=0.1
    )
    
    # Initial conditions
    n_grid = nx * ny
    u0 = zeros(2 + n_grid + 4)
    u0[1] = -65.0
    u0[2] = 0.0
    
    # Time span
    tspan = (0.0, 1000.0)
    
    # Create and solve ODE problem
    prob = ODEProblem((du, u, p, t) -> dopaminergic_autoreceptor_grid!(du, u, p, t, grid_params, :fdm),
                     u0, tspan, p)
    sol = solve(prob, Tsit5(), saveat=10.0)
    
    # Plot results using GLMakie
    fig = Figure(size=(1200, 400))
    times = [100, 500, 1000]
    
    for (i, t) in enumerate(times)
        idx = findfirst(x -> x >= t, sol.t)
        DA_grid = reshape(sol.u[idx][3:end-4], nx, ny)
        
        ax = Axis(fig[1, i], title="t = $t")
        heatmap!(ax, DA_grid, colormap=:viridis)
        Colorbar(fig[1, i+1], limits=(0, maximum(DA_grid)))
    end
    
    display(fig)
end

"""
    compare_methods()
Compare the performance of different discretization methods.
"""
function compare_methods()
    # Grid parameters
    nx, ny = 50, 50
    dx, dy = 1.0, 1.0
    D = 0.1
    release_sites = [(25,25)]
    
    grid_params = GridParameters(nx, ny, dx, dy, D, release_sites)
    
    # Model parameters
    p = (
        Cm=1.0, gL=0.1, EL=-65.0, gCa=1.0, ECa=120.0,
        Vhalf=-20.0, kV=10.0, τCa=100.0, αCa=0.1,
        krel=1.0, kclear=0.1, kon=1.0, koff=0.1,
        kG=1.0, kGdeg=0.1, G50=1.0, nGi=2.0,
        kACbasal=1.0, kcAMPdeg=0.1, kPKA=1.0, kPKAdeg=0.1
    )
    
    # Initial conditions
    n_grid = nx * ny
    u0 = zeros(2 + n_grid + 4)
    u0[1] = -65.0
    u0[2] = 0.0
    
    # Time span
    tspan = (0.0, 100.0)  # Shorter time span for comparison
    
    # Compare methods
    methods = [:fdm, :fvm, :spectral]
    times = Dict{Symbol,Float64}()
    
    for method in methods
        prob = ODEProblem((du, u, p, t) -> dopaminergic_autoreceptor_grid!(du, u, p, t, grid_params, method),
                         u0, tspan, p)
        t = @elapsed solve(prob, Tsit5())
        times[method] = t
    end
    
    # Print results
    println("Method comparison:")
    for (method, time) in times
        println("$method: $time seconds")
    end
end

# Run tests
println("Testing single release site...")
test_single_release_site()

println("\nTesting multiple release sites...")
test_multiple_release_sites()

println("\nComparing methods...")
compare_methods() 
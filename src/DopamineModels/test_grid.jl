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

# Run the model with the dopamine grid Initial conditions
n_grid = nx * ny
u0 = zeros(2 + n_grid + 4)  # V, Ca, DA_grid, P, Gi, cAMP, PKA
u0[1] = -65.0  # Initial voltage
u0[2] = 0.0    # Initial calcium

# Time span
tspan = (0.0, 1000.0)

# Create the ODE problem
prob = ODEProblem(
    (du, u, p, t) -> dopaminergic_autoreceptor_grid!(du, u, p, t, grid_params; method = :fdm),
                    u0, tspan, p)

# Solve
sol = solve(prob, Tsit5(), saveat=10.0)

# Create 3D array of solutions
DA_grid = cat(map(t -> reshape(sol(t)[3:end-4], nx, ny), sol.t)...; dims=3)

#%% Create the animation
fig = Figure(size=(1200, 800))

# Left panel for heatmap
ax1 = Axis(fig[1, 1], title="Dopamine Concentration Evolution", aspect=1.0)
# Calculate the maximum value across all frames for consistent color limits
max_val = maximum(DA_grid)
hm = heatmap!(ax1, DA_grid[:,:,1], colormap=:viridis, colorrange=(0, max_val))
Colorbar(fig[1, 2], hm)

# Add release site markers
for (i, j) in release_sites
    scatter!(ax1, [i], [j], color=:red, markersize=15, label="Release Site")
end

# Set axis limits and labels
xlims!(ax1, 0.5, nx+0.5)
ylims!(ax1, 0.5, ny+0.5)
axislegend(ax1)

# Right panel for time traces
ax2 = Axis(fig[1, 3], title="Solution Values", xlabel="Time", ylabel="Value")
# Create lines for each variable
lines!(ax2, sol.t, [u[1] for u in sol.u], label="V", color=:blue)
lines!(ax2, sol.t, [u[2] for u in sol.u], label="Ca", color=:red)
lines!(ax2, sol.t, [mean(reshape(u[3:end-4], nx, ny)) for u in sol.u], label="DA_avg", color=:green)
lines!(ax2, sol.t, [u[end-3] for u in sol.u], label="P", color=:purple)
lines!(ax2, sol.t, [u[end-2] for u in sol.u], label="Gi", color=:orange)
lines!(ax2, sol.t, [u[end-1] for u in sol.u], label="cAMP", color=:brown)
lines!(ax2, sol.t, [u[end] for u in sol.u], label="PKA", color=:pink)

# Add a vertical line to show current time
vline = lines!(ax2, [0], [0], color=:black, linewidth=2)
axislegend(ax2)
fig
#%%
# Create animation
framerate = 30
record(fig, "src/DopamineModels/TestImages/dopamine_evolution.mp4", 1:size(DA_grid, 3); framerate=framerate) do frame
    println(frame)
    hm[1] = DA_grid[:,:,frame]  # Update the heatmap data
    ax1.title = "Dopamine Concentration at t = $(round(sol.t[frame], digits=1))"
    
    # Update the vertical line position
    vline[1] = [sol.t[frame], sol.t[frame]]
    vline[2] = [minimum(ax2.ylims[]), maximum(ax2.ylims[])]
end
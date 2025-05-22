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

# -Test the grid-based model with a single release site in the center.

# Grid parameters
nx, ny = 50, 50
dx, dy = 1.0, 1.0
D = 0.01  # Diffusion coefficient
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
sol = solve(prob, Tsit5(), saveat=1.0)

# Create 3D array of solutions
DA_grid = cat(map(t -> reshape(sol(t)[3:end-4], nx, ny), sol.t)...; dims=3)

# Create the animation
fig = Figure(size=(1200, 1000))

# Left panel for heatmap
ax1 = Axis(fig[1:3, 1], title="Dopamine Concentration Evolution", aspect=1.0)
# Calculate the maximum value across all frames for consistent color limits
max_val = maximum(DA_grid)
hm = heatmap!(ax1, DA_grid[:,:,1], colormap=:viridis, colorrange=(0, max_val))
Colorbar(fig[1:3, 2], hm)

# Add release site markers
for (i, j) in release_sites
    scatter!(ax1, [i], [j], color=:red, markersize=15, label="Release Site")
end

# Set axis limits and labels
xlims!(ax1, 0.5, nx+0.5)
ylims!(ax1, 0.5, ny+0.5)
axislegend(ax1)

# Right panel for time traces - one subplot per variable
ax_v = Axis(fig[1, 3], title="Voltage", xlabel="Time", ylabel="V (mV)")
ax_ca = Axis(fig[2, 3], title="Calcium", xlabel="Time", ylabel="Ca (μM)")
ax_da = Axis(fig[3, 3], title="Average Dopamine", xlabel="Time", ylabel="DA (μM)")

# Create lines for each variable
lines!(ax_v, sol.t, [u[1] for u in sol.u], color=:blue)
lines!(ax_ca, sol.t, [u[2] for u in sol.u], color=:red)
lines!(ax_da, sol.t, [mean(reshape(u[3:end-4], nx, ny)) for u in sol.u], color=:green)

# Add vertical lines to show current time
vline_v = vlines!(ax_v, [0], color=:black, linewidth=2)
vline_ca = vlines!(ax_ca, [0], color=:black, linewidth=2)
vline_da = vlines!(ax_da, [0], color=:black, linewidth=2)

# Link x-axes of all time trace plots
linkxaxes!(ax_v, ax_ca, ax_da)

# Create animation
framerate = 30
record(fig, "src/DopamineModels/TestImages/dopamine_evolution.mp4", 1:size(DA_grid, 3); framerate=framerate) do frame
    println(frame)
    hm[1] = DA_grid[:,:,frame]  # Update the heatmap data
    ax1.title = "Dopamine Concentration at t = $(round(sol.t[frame], digits=1))"
    
    # Update the vertical line positions
    current_time = sol.t[frame]
    vline_v[1] = [current_time]
    vline_ca[1] = [current_time]
    vline_da[1] = [current_time]
end

fig
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

# ---------- Create the parameters for the dopamine grid ----------
release_sites = [(5, 5)]  # Center of the grid
grid_params = GridParameters(
    settings.nx, settings.ny, 
    settings.dx, settings.dy, 
    settings.D, release_sites
)
grid_params.fdm_operator

p_grid = (p.krel, p.kclear)

#Create the initial conditions for the dopamine grid
u0 = zeros(settings.nx*settings.ny)
# Calculate middle point indices
mid_x = div(settings.nx, 2) + 1  # Integer division by 2, then add 1
mid_y = div(settings.ny, 2) + 1
# Convert 2D index to 1D index
mid_idx = (mid_y - 1) * settings.nx + mid_x
u0[mid_idx] = 1.0

tspan = (0.0, 100.0)

prob = ODEProblem(
    (du, u, p, t) -> update_dopamine_grid!(du, u, p, t, grid_params),
    u0, tspan, p_grid
)

sol = solve(prob, Tsit5(), saveat=1.0)  # Save more frequently for smoother animation

# Create 3D array of solutions
DA_grid = cat(map(t -> reshape(sol(t), settings.nx, settings.ny), sol.t)...; dims=3)

# Create animation
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
xlims!(ax, 0.5, settings.nx+0.5)
ylims!(ax, 0.5, settings.ny+0.5)
axislegend(ax)

# Create animation
framerate = 5
record(fig, "src/DopamineModels/TestImages/grid_diffusion_evolution.mp4", 1:size(DA_grid, 3); framerate=framerate) do frame
    println(frame)
    hm[1] = DA_grid[:,:,frame]  # Update the heatmap data
    ax.title = "Dopamine Concentration at t = $(round(sol.t[frame], digits=1))"
end
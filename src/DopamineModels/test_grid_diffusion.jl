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
u0[mid_idx] = 10.0

tspan = (0.0, 100.0)

prob = ODEProblem(
    (du, u, p, t) -> update_dopamine_grid!(du, u, p, t, grid_params),
    u0, tspan, p_grid
)

sol = solve(prob, Tsit5(), saveat=1.0)  # Save more frequently for smoother animation

# Create 3D array of solutions
DA_grid = cat(map(t -> reshape(sol(t), settings.nx, settings.ny), sol.t)...; dims=3)

# Create animation
fig = Figure(size=(1200, 1000))

# Left panel for heatmap
ax1 = Axis(fig[1:2, 1], title="Dopamine Concentration Evolution", aspect=1.0)
# Calculate the maximum value across all frames for consistent color limits
max_val = maximum(DA_grid)
hm = heatmap!(ax1, DA_grid[:,:,1], colormap=:viridis, colorrange=(0, max_val))
Colorbar(fig[1:2, 0], hm)

# Add release site markers
for (i, j) in release_sites
    scatter!(ax1, [i], [j], color=:red, markersize=15, label="Release Site")
end

# Set axis limits and labels
xlims!(ax1, 0.5, settings.nx+0.5)
ylims!(ax1, 0.5, settings.ny+0.5)
axislegend(ax1)

# Add subplots for X and Y line averages
ax_x = Axis(fig[3, 1], title="Average DA along X", xlabel="X position", ylabel="DA (μM)")
ax_y = Axis(fig[1:2, 2], title="Average DA along Y", xlabel="Y position", ylabel="DA (μM)")

# Create lines for X and Y averages
x_line = lines!(ax_x, 1:settings.nx, zeros(settings.nx), color=:blue)
y_line = lines!(ax_y, zeros(settings.ny), 1:settings.ny,color=:red)

# Create animation
framerate = 30
record(fig, "src/DopamineModels/TestImages/grid_diffusion.mp4", 1:size(DA_grid, 3); framerate=framerate) do frame
    println(frame)
    current_grid = DA_grid[:,:,frame]
    hm[1] = current_grid  # Update the heatmap data
    ax1.title = "Dopamine Concentration at t = $(round(sol.t[frame], digits=1))"
    
    # Update X and Y averages
    x_line[2] = mean(current_grid, dims=1)[:]  # Average along Y for each X
    y_line[1] = mean(current_grid, dims=2)[:]  # Average along X for each Y
end
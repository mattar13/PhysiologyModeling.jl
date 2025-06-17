using DifferentialEquations
using SparseArrays
using LinearAlgebra
using Plots
using ProgressMeter
using Plots

# --------------------------------------------------
#%% 3) Initialize and solve with progress bar
# --------------------------------------------------
include("../src/MitoModel/auxillary_functions.jl")
include("../src/MitoModel/parameters.jl")
include("../src/MitoModel/models.jl")

prob = ODEProblem(
    (du, u, p, t) -> pde_system!(du, u, p, t; N = N), 
    u0, tspan, params
)

# use a solver with constant Jacobian and show progress
sol = solve(prob, Rodas5(autodiff=false), progress=true)

# --------------------------------------------------
#%% 4) Visualization of initial & final states
# --------------------------------------------------
# Function to add mitochondria outline and contour to a plot
function add_mito_visualization!(p, params::Params2D)
    # Add outline
    #plot!(p, x, y, linecolor=:red, linewidth=2, label="")
    # Add contour of mitochondrial region
    mito_mask_2d = reshape(params.mito_mask, params.ny, params.nx)
    contour!(p, params.xrng, params.yrng, mito_mask_2d, 
        fillalpha=1.0,  # Very low alpha for subtle effect
        levels=[0.0, 0.01],  # Only show the boundary between 0 and 1
        color=:red,
        label="")
end

# Initial conditions
p1 = heatmap(params.xrng, params.yrng, reshape(u0[1:N], params.ny, params.nx);
    title="Initial Ca²⁺ Distribution", aspect_ratio=1,
    xlabel="x", ylabel="y", cbar_title="[Ca]")
add_mito_visualization!(p1, params)

p2 = heatmap(params.xrng, params.yrng, reshape(u0[N+1:2N], params.ny, params.nx);
    title="Initial ATP Distribution", aspect_ratio=1,
    xlabel="x", ylabel="y", cbar_title="[ATP]")
add_mito_visualization!(p2, params)

p3 = heatmap(params.xrng, params.yrng, reshape(u0[2N+1:3N], params.ny, params.nx);
    title="Initial ADP Distribution", aspect_ratio=1,
    xlabel="x", ylabel="y", cbar_title="[ADP]")
add_mito_visualization!(p3, params)

p4 = heatmap(params.xrng, params.yrng, reshape(u0[3N+1:4N], params.ny, params.nx);
    title="Initial AMP Distribution", aspect_ratio=1,
    xlabel="x", ylabel="y", cbar_title="[AMP]")
add_mito_visualization!(p4, params)

p5 = heatmap(params.xrng, params.yrng, reshape(u0[4N+1:5N], params.ny, params.nx);
    title="Initial Adenosine Distribution", aspect_ratio=1,
    xlabel="x", ylabel="y", cbar_title="[Adenosine]")
add_mito_visualization!(p5, params)

p6 = heatmap(params.xrng, params.yrng, reshape(u0[5N+1:6N], params.ny, params.nx);
    title="Initial Phosphate Distribution", aspect_ratio=1,
    xlabel="x", ylabel="y", cbar_title="[P]")
add_mito_visualization!(p6, params)

plot(p1, p2, p3, p4, p5, p6, layout=(2,3), size=(1200,800))

# Final states
p7 = heatmap(params.xrng, params.yrng, reshape(sol.u[end][1:N], params.ny, params.nx);
    title="Final Ca²⁺ (t=$(round(sol.t[end],digits=1)) ms)", aspect_ratio=1,
    xlabel="x", ylabel="y", cbar_title="[Ca]")
add_mito_visualization!(p7, params)

p8 = heatmap(params.xrng, params.yrng, reshape(sol.u[end][N+1:2N], params.ny, params.nx);
    title="Final ATP (t=$(round(sol.t[end],digits=1)) ms)", aspect_ratio=1,
    xlabel="x", ylabel="y", cbar_title="[ATP]")
add_mito_visualization!(p8, params)

p9 = heatmap(params.xrng, params.yrng, reshape(sol.u[end][2N+1:3N], params.ny, params.nx);
    title="Final ADP (t=$(round(sol.t[end],digits=1)) ms)", aspect_ratio=1,
    xlabel="x", ylabel="y", cbar_title="[ADP]")
add_mito_visualization!(p9, params)

p10 = heatmap(params.xrng, params.yrng, reshape(sol.u[end][3N+1:4N], params.ny, params.nx);
    title="Final AMP (t=$(round(sol.t[end],digits=1)) ms)", aspect_ratio=1,
    xlabel="x", ylabel="y", cbar_title="[AMP]")
add_mito_visualization!(p10, params)

p11 = heatmap(params.xrng, params.yrng, reshape(sol.u[end][4N+1:5N], params.ny, params.nx);
    title="Final Adenosine (t=$(round(sol.t[end],digits=1)) ms)", aspect_ratio=1,
    xlabel="x", ylabel="y", cbar_title="[Adenosine]")
add_mito_visualization!(p11, params)

p12 = heatmap(params.xrng, params.yrng, reshape(sol.u[end][5N+1:6N], params.ny, params.nx);
    title="Final Phosphate (t=$(round(sol.t[end],digits=1)) ms)", aspect_ratio=1,
    xlabel="x", ylabel="y", cbar_title="[P]")
add_mito_visualization!(p12, params)

plot(p7, p8, p9, p10, p11, p12, layout=(2,3), size=(1200,800))

# --------------------------------------------------
#%% 5) Create animation of the diffusion process
# --------------------------------------------------
# Create animation
anim = @animate for i in 1:length(sol.t)
    # Get current state
    ca = reshape(sol.u[i][1:N], ny, nx)
    atp = reshape(sol.u[i][N+1:2N], ny, nx)
    adp = reshape(sol.u[i][2N+1:3N], ny, nx)
    amp = reshape(sol.u[i][3N+1:4N], ny, nx)
    ado = reshape(sol.u[i][4N+1:5N], ny, nx)
    p = reshape(sol.u[i][5N+1:6N], ny, nx)
    
    # Create plots
    p1 = heatmap(xrng, yrng, ca, 
        title="Ca²⁺ (t=$(round(sol.t[i],digits=1)) ms)",
        aspect_ratio=1,
        xlabel="x", ylabel="y", 
        cbar_title="[Ca]",
        clims=(0, 0.1))
    add_mito_visualization!(p1, params)
    
    p2 = heatmap(xrng, yrng, atp,
        title="ATP (t=$(round(sol.t[i],digits=1)) ms)",
        aspect_ratio=1,
        xlabel="x", ylabel="y",
        cbar_title="[ATP]",
        clims=(0, 0.1))
    add_mito_visualization!(p2, params)
    
    p3 = heatmap(xrng, yrng, adp,
        title="ADP (t=$(round(sol.t[i],digits=1)) ms)",
        aspect_ratio=1,
        xlabel="x", ylabel="y",
        cbar_title="[ADP]",
        clims=(0, 0.001))
    add_mito_visualization!(p3, params)
    
    p4 = heatmap(xrng, yrng, amp,
        title="AMP (t=$(round(sol.t[i],digits=1)) ms)",
        aspect_ratio=1,
        xlabel="x", ylabel="y",
        cbar_title="[AMP]",
        clims=(0, 0.001))
    add_mito_visualization!(p4, params)
    
    p5 = heatmap(xrng, yrng, ado,
        title="Adenosine (t=$(round(sol.t[i],digits=1)) ms)",
        aspect_ratio=1,
        xlabel="x", ylabel="y",
        cbar_title="[Adenosine]",
        clims=(0, 0.001))
    add_mito_visualization!(p5, params)
    
    p6 = heatmap(xrng, yrng, p,
        title="Phosphate (t=$(round(sol.t[i],digits=1)) ms)",
        aspect_ratio=1,
        xlabel="x", ylabel="y",
        cbar_title="[P]",
        clims=(0, 0.1))
    add_mito_visualization!(p6, params)
    
    # Combine plots
    plot(p1, p2, p3, p4, p5, p6, layout=(2,3), size=(1200,800))
end

# Save animation
gif(anim, "diffusion_animation.gif", fps=15)

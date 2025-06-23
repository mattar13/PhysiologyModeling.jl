using DifferentialEquations
using SparseArrays
using LinearAlgebra
using Plots
using ProgressMeter
using Plots

# --------------------------------------------------
#%% 3) Initialize and solve with progress bar
# --------------------------------------------------
println("Initializing...")
include("../src/MitoModel/auxillary_functions.jl")
include("../src/MitoModel/parameters.jl")
include("../src/MitoModel/models.jl")

println("Solving...")
# Create ODE problem
prob = ODEProblem(
    (du, u, p, t) -> pde_system!(du, u, p, t; N = N), 
    u0, tspan, params
)

# use a solver with constant Jacobian and show progress
# Add tstops to ensure we capture stimulus start and end
sol = solve(prob, Rodas5(autodiff=false), 
    tstops=[stim_start, stim_end],
    progress=true
)


# --------------------------------------------------
#%% 5) Create animation of the diffusion process
# --------------------------------------------------
# Function to add mitochondria outline and contour to a plot
function add_mito_visualization!(p, params::Params2D)
    # Add contour of mitochondrial region (red)
    mito_mask_2d = reshape(params.mito_mask .+ params.cyto_mask, params.ny, params.nx)
    # heatmap!(p, params.xrng, params.yrng, mito_mask_2d./10.0, 
    #     alpha = 0.01, color = :greens # This will clip values below 0.5
    # )
end

nt = 100
t_rng = LinRange(0.0, sol.t[end], nt)
sol_arr = Array(sol(t_rng))

ca_t = reshape(sol_arr[1:N, :], params.ny, params.nx, nt)
atp_t = reshape(sol_arr[N+1:2N, :], params.ny, params.nx, nt)
adp_t = reshape(sol_arr[2N+1:3N, :], params.ny, params.nx, nt)
amp_t = reshape(sol_arr[3N+1:4N, :], params.ny, params.nx, nt)
ado_t = reshape(sol_arr[4N+1:5N, :], params.ny, params.nx, nt)
p_t = reshape(sol_arr[5N+1:6N, :], params.ny, params.nx, nt)
Vm_t = reshape(sol_arr[6N+1:7N, :], params.ny, params.nx, nt)


# Create animation
anim = @animate for i in eachindex(t_rng)
    t = t_rng[i]
    println("t = $t")
    # Get current state
    Vm_i = Vm_t[:, :, i]
    atp_i = atp_t[:, :, i]
    adp_i = adp_t[:, :, i]
    amp_i = amp_t[:, :, i]
    ado_i = ado_t[:, :, i]
    p_i = p_t[:, :, i]
    
    # Create plots
    p1 = heatmap(params.xrng, params.yrng, Vm_i, 
        title="Vm (t=$(round(t ,digits=1)) ms)",
        aspect_ratio=1,
        xlabel="x", ylabel="y", 
        cbar_title="[Vm]",
        clims=(-100, 0)
    )
    add_mito_visualization!(p1, params)
    
    p2 = heatmap(params.xrng, params.yrng, atp_i,
        title="ATP (t=$(round(t ,digits=1)) ms)",
        aspect_ratio=1,
        xlabel="x", ylabel="y",
        cbar_title="[ATP]",
        clims=(0, maximum(atp_t))
    )
    add_mito_visualization!(p2, params)
    
    p3 = heatmap(params.xrng, params.yrng, adp_i,
        title="ADP (t=$(round(t ,digits=1)) ms)",
        aspect_ratio=1,
        xlabel="x", ylabel="y",
        cbar_title="[ADP]",
        clims=(0, maximum(adp_t))
    )
    add_mito_visualization!(p3, params)
    
    p4 = heatmap(params.xrng, params.yrng, amp_i,
        title="AMP (t=$(round(t ,digits=1)) ms)",
        aspect_ratio=1,
        xlabel="x", ylabel="y",
        cbar_title="[AMP]",
        clims=(0, maximum(amp_t))
    )
    add_mito_visualization!(p4, params)
    
    p5 = heatmap(params.xrng, params.yrng, ado_i,
        title="Adenosine (t=$(round(t ,digits=1)) ms)",
        aspect_ratio=1,
        xlabel="x", ylabel="y",
        cbar_title="[Adenosine]",
        clims=(0, maximum(ado_t))
    )
    add_mito_visualization!(p5, params)
    
    p6 = heatmap(params.xrng, params.yrng, p_i ,
        title="Phosphate (t=$(round(t ,digits=1)) ms)",
        aspect_ratio=1,
        xlabel="x", ylabel="y",
        cbar_title="[P]",
        clims=(0, maximum(p_t))
    )
    add_mito_visualization!(p6, params)
    
    # Combine plots
    plot(p1, p2, p3, p4, p5, p6, layout=(2,3), size=(1200,800))
end

# Save animation
gif(anim, "diffusion_animation.gif", fps=15)


p1 = plot(Vm_t[5,5,:])
p2 = plot(atp_t[5,5,:])
plot(p1, p2, layout=(2,1), size=(1200,800))
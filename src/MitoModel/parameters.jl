# --------------------------------------------------
#%% 1) Define parameters
# --------------------------------------------------

# Package parameters into struct
struct Params2D
    #Diffusion coefficients
    nx::Int
    ny::Int
    xrng::Vector{Float64}
    yrng::Vector{Float64}
    D_ca::Float64
    D_atp::Float64
    D_adp::Float64
    D_amp::Float64
    D_ado::Float64
    D_p::Float64    # Diffusion coefficient for phosphate
    #Reaction rates
    k_atp_adp::Float64
    k_adp_amp::Float64
    k_amp_ado::Float64
    k_adk::Float64
    k_ak::Float64  # Rate constant for mitochondrial adenylate kinase reaction
    k_atp_synth::Float64  # Rate constant for ATP synthesis from ADP + P
    k_atp_prod::Float64
    #Laplacian matrix
    L::SparseMatrixCSC{Float64,Int}
    #Mitochondria mask
    mito_mask::Vector{Bool}
end

# --------------------------------------------------
#%% 1) Define grid parameters, and precompute Laplacian
# --------------------------------------------------
nx, ny = 100, 100            # Grid resolution
xrng = LinRange(0, 1, nx)
yrng = LinRange(0, 1, ny)
dx = xrng[2] - xrng[1]
dy = yrng[2] - yrng[1]
N = nx*ny
# Precompute Laplacian matrix
L = build_laplacian(nx, ny, dx, dy)

# --------------------------------------------------
#%% 2) Define parameters
# --------------------------------------------------
D_ca = 0.1                  # Diffusion coefficient for Ca²⁺ (µm²/ms)
D_atp = 0.3                 # Diffusion coefficient for ATP (µm²/ms)
D_adp = 0.3                 # Diffusion coefficient for ADP (µm²/ms)
D_amp = 0.3                 # Diffusion coefficient for AMP (µm²/ms)
D_ado = 0.3                 # Diffusion coefficient for Adenosine (µm²/ms)
D_p = 0.3                   # Diffusion coefficient for Phosphate (µm²/ms)
k_atp_adp = 0.01           # Rate constant for ATP → ADP conversion (ms⁻¹)
k_adp_amp = 0.008          # Rate constant for ADP → AMP conversion (ms⁻¹)
k_amp_ado = 0.005          # Rate constant for AMP → Adenosine conversion (ms⁻¹)
k_adk = 0.02               # Rate constant for ADK reaction: Adenosine + ATP → AMP + ADP (ms⁻¹)
k_ak = 0.015               # Rate constant for mitochondrial AK reaction: 2ADP → ATP + AMP (ms⁻¹)
k_atp_synth = 2.0         # Rate constant for ATP synthesis: ADP + P → ATP (ms⁻¹)
k_atp_prod = 0.05          # Rate constant for ATP production in mitochondria (ms⁻¹)
tspan = (0.0, 50.0)        # Time span
N = nx * ny                 # total grid points


# Define mitochondria location and size
mito_center = (0.5, 0.5)    # Center of mitochondria (normalized coordinates)
mito_radius = 0.10          # Radius of mitochondria (normalized units)

# Create mitochondria mask
mito_mask = create_mito_mask(nx, ny, mito_center, mito_radius, dx = dx, dy = dy)


# initial condition: uniform concentration for all species
u0 = Float64[]
ca0 = rand(Float64, N)*0.1
atp0 = ones(Float64, N)*3.0#; atp0[25] = 1.0
adp0 = ones(Float64, N)*0.0
amp0 = ones(Float64, N)*0.0
ado0 = ones(Float64, N)*0.0
p0 = ones(Float64, N)*0.0  # Initial phosphate concentration
push!(u0, ca0...)
push!(u0, atp0...)
push!(u0, adp0...)
push!(u0, amp0...)
push!(u0, ado0...)
push!(u0, p0...)

# pack parameters and define ODE problem
params = Params2D(nx, ny, xrng, yrng, 
    D_ca, D_atp, D_adp, D_amp, D_ado, D_p,
    k_atp_adp, k_adp_amp, k_amp_ado, k_adk, k_ak, k_atp_synth, k_atp_prod,
    L, mito_mask
)
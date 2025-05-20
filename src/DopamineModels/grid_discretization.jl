using SparseArrays
using LinearAlgebra

"""
    GridParameters
Parameters for the 2D grid discretization
"""
struct GridParameters
    nx::Int          # Number of points in x direction
    ny::Int          # Number of points in y direction
    dx::Float64      # Grid spacing in x
    dy::Float64      # Grid spacing in y
    D::Float64       # Diffusion coefficient
    release_sites::Vector{Tuple{Int,Int}}  # Coordinates of release sites
end

"""
    create_fdm_operator(params::GridParameters)
Creates the finite difference method operator for 2D diffusion
"""
function create_fdm_operator(params::GridParameters)
    nx, ny = params.nx, params.ny
    dx, dy = params.dx, params.dy
    D = params.D
    
    # Create sparse matrix for the Laplacian
    n = nx * ny
    I = Int[]
    J = Int[]
    V = Float64[]
    
    for j in 1:ny, i in 1:nx
        idx = (j-1)*nx + i
        # Center point
        push!(I, idx)
        push!(J, idx)
        push!(V, -2D*(1/dx^2 + 1/dy^2))
        
        # Neighbors
        if i > 1
            push!(I, idx)
            push!(J, idx-1)
            push!(V, D/dx^2)
        end
        if i < nx
            push!(I, idx)
            push!(J, idx+1)
            push!(V, D/dx^2)
        end
        if j > 1
            push!(I, idx)
            push!(J, idx-nx)
            push!(V, D/dy^2)
        end
        if j < ny
            push!(I, idx)
            push!(J, idx+nx)
            push!(V, D/dy^2)
        end
    end
    
    return sparse(I, J, V, n, n)
end

"""
    create_fvm_operator(params::GridParameters)
Creates the finite volume method operator for 2D diffusion
"""
function create_fvm_operator(params::GridParameters)
    nx, ny = params.nx, params.ny
    dx, dy = params.dx, params.dy
    D = params.D
    
    n = nx * ny
    I = Int[]
    J = Int[]
    V = Float64[]
    
    for j in 1:ny, i in 1:nx
        idx = (j-1)*nx + i
        # Center point
        push!(I, idx)
        push!(J, idx)
        push!(V, -2D*(1/dx^2 + 1/dy^2))
        
        # Face fluxes
        if i > 1
            push!(I, idx)
            push!(J, idx-1)
            push!(V, D/dx^2)
        end
        if i < nx
            push!(I, idx)
            push!(J, idx+1)
            push!(V, D/dx^2)
        end
        if j > 1
            push!(I, idx)
            push!(J, idx-nx)
            push!(V, D/dy^2)
        end
        if j < ny
            push!(I, idx)
            push!(J, idx+nx)
            push!(V, D/dy^2)
        end
    end
    
    return sparse(I, J, V, n, n)
end

"""
    create_spectral_operator(params::GridParameters)
Creates the spectral method operator for 2D diffusion
"""
function create_spectral_operator(params::GridParameters)
    nx, ny = params.nx, params.ny
    dx, dy = params.dx, params.dy
    D = params.D
    
    # Create wave numbers
    kx = 2π * fftfreq(nx, 1/dx)
    ky = 2π * fftfreq(ny, 1/dy)
    
    # Create the spectral operator
    K = zeros(Complex{Float64}, nx, ny)
    for j in 1:ny, i in 1:nx
        K[i,j] = -D * (kx[i]^2 + ky[j]^2)
    end
    
    return K
end

"""
    update_dopamine_grid!(du, u, p, t, grid_params, method)
Updates the dopamine concentration grid using the specified method
"""
function update_dopamine_grid!(du, u, p, t, grid_params, method)
    krel, kclear = p[3], p[4]  # Get release and clearance parameters
    nx, ny = grid_params.nx, grid_params.ny
    
    # Reshape the dopamine grid
    DA_grid = reshape(u[3:end-4], nx, ny)
    
    # Add release at specific sites
    for (i, j) in grid_params.release_sites
        if 1 ≤ i ≤ nx && 1 ≤ j ≤ ny
            DA_grid[i,j] += krel * u[2]  # Add release based on calcium
        end
    end
    
    # Apply diffusion operator based on method
    if method == :fdm
        L = create_fdm_operator(grid_params)
        du[3:end-4] = L * vec(DA_grid) - kclear * vec(DA_grid)
    elseif method == :fvm
        L = create_fvm_operator(grid_params)
        du[3:end-4] = L * vec(DA_grid) - kclear * vec(DA_grid)
    elseif method == :spectral
        K = create_spectral_operator(grid_params)
        DA_hat = fft(DA_grid)
        DA_hat .*= exp.(K * dt)
        du[3:end-4] = real(ifft(DA_hat)) - kclear * vec(DA_grid)
    end
end

"""
    benchmark_methods(grid_params, tspan, u0, p)
Benchmarks the different discretization methods
"""
function benchmark_methods(grid_params, tspan, u0, p)
    methods = [:fdm, :fvm, :spectral]
    times = Dict{Symbol,Float64}()
    
    for method in methods
        t = @elapsed begin
            # Solve the system using the specified method
            # (Implementation depends on your ODE solver)
        end
        times[method] = t
    end
    
    return times
end

using SparseArrays
using LinearAlgebra
using FFTW

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
    fdm_operator::SparseMatrixCSC{Float64,Int}  # FDM operator
    fvm_operator::SparseMatrixCSC{Float64,Int}  # FVM operator
    spectral_operator::Matrix{Complex{Float64}}  # Spectral operator
end

"""
    GridParameters(nx, ny, dx, dy, D, release_sites)
Constructor for GridParameters that creates and stores the diffusion operators.
"""
function GridParameters(nx::Int, ny::Int, dx::Float64, dy::Float64, D::Float64, release_sites::Vector{Tuple{Int,Int}})
    fdm_op = create_fdm_operator(nx, ny, dx, dy, D)
    fvm_op = create_fvm_operator(nx, ny, dx, dy, D)
    spec_op = create_spectral_operator(nx, ny, dx, dy, D)
    
    return GridParameters(nx, ny, dx, dy, D, release_sites, fdm_op, fvm_op, spec_op)
end

"""
    create_fdm_operator(nx, ny, dx, dy, D)
Creates the finite difference method operator for 2D diffusion
"""
function create_fdm_operator(nx::Int, ny::Int, dx::Float64, dy::Float64, D::Float64)
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
    create_fvm_operator(nx, ny, dx, dy, D)
Creates the finite volume method operator for 2D diffusion
"""
function create_fvm_operator(nx::Int, ny::Int, dx::Float64, dy::Float64, D::Float64)
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
    create_spectral_operator(nx, ny, dx, dy, D)
Creates the spectral method operator for 2D diffusion
"""
function create_spectral_operator(nx::Int, ny::Int, dx::Float64, dy::Float64, D::Float64)
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
    update_dopamine_grid!(dDA_grid, DA_grid, p, t, grid_params; method=:fdm)
Updates the dopamine concentration grid using the specified method.
Inputs:
- dDA_grid: Vector for storing derivatives
- DA_grid: Current dopamine concentration vector
- p: Tuple of (krel, kclear)
- t: Current time
- grid_params: Grid parameters
- method: Discretization method (:fdm, :fvm, or :spectral)
"""
function update_dopamine_grid!(dDA_grid::Vector{T}, DA_grid::Vector{T}, p, t, grid_params; method=:fdm) where T <: Real
    krel, kclear = p  # Get release and clearance parameters
    nx, ny = grid_params.nx, grid_params.ny
    
    # Add release at specific sites
    # for (i, j) in grid_params.release_sites
    #     if 1 ≤ i ≤ nx && 1 ≤ j ≤ ny
    #         idx = (j-1)*nx + i
    #         DA_grid[idx] += krel * t  # Add release based on calcium (t is actually Ca here)
    #     end
    # end
    
    # Apply diffusion operator based on method
    if method == :fdm
        dDA_grid .= grid_params.fdm_operator * DA_grid - kclear * DA_grid
    elseif method == :fvm
        dDA_grid .= grid_params.fvm_operator * DA_grid - kclear * DA_grid
    elseif method == :spectral
        K = grid_params.spectral_operator
        DA_hat = fft(reshape(DA_grid, nx, ny))
        DA_hat .*= exp.(K * dt)
        dDA_grid .= vec(real(ifft(DA_hat)))# - kclear * DA_grid
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

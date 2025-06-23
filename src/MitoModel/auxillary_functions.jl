α_m(V) = 0.1 * (25.0 - V) / (exp((25.0 - V)/10.0) - 1.0)
β_m(V) = 4.0 * exp(-V/18.0)
α_h(V) = 0.07 * exp(-V/20.0)
β_h(V) = 1.0 / (exp((30.0 - V)/10.0) + 1.0)
α_n(V) = 0.01 * (10.0 - V) / (exp((10.0 - V)/10.0) - 1.0)
β_n(V) = 0.125 * exp(-V/80.0)


"""
Build a sparse 2D Laplacian with Neumann BC via reflection:
- Diagonal: -2/dx² -2/dy²
- Off-diagonals: +1/dx² (left/right), +1/dy² (up/down)
"""
function build_laplacian(nx::Int, ny::Int, dx::Float64, dy::Float64)
    N = nx * ny
    invdx2 = 1.0 / dx^2
    invdy2 = 1.0 / dy^2
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    for j in 1:ny, i in 1:nx
        k = (j-1)*nx + i
        # center
        push!(rows, k); push!(cols, k); push!(vals, -2*invdx2 - 2*invdy2)
        # left neighbor (Neumann BC reflection)
        il = i > 1 ? i-1 : 2-i
        kL = (j-1)*nx + il
        push!(rows, k); push!(cols, kL); push!(vals, invdx2)
        # right neighbor
        ir = i < nx ? i+1 : 2*nx - i
        kR = (j-1)*nx + ir
        push!(rows, k); push!(cols, kR); push!(vals, invdx2)
        # down neighbor
        jd = j > 1 ? j-1 : 2-j
        kD = (jd-1)*nx + i
        push!(rows, k); push!(cols, kD); push!(vals, invdy2)
        # up neighbor
        ju = j < ny ? j+1 : 2*ny - j
        kU = (ju-1)*nx + i
        push!(rows, k); push!(cols, kU); push!(vals, invdy2)
    end
    return sparse(rows, cols, vals, N, N)
end

"""
Build a sparse 2D Laplacian with mask-based boundary conditions:
- Diagonal: -2/dx² -2/dy²
- Off-diagonals: +1/dx² (left/right), +1/dy² (up/down)
- Zero flux at masked boundaries
"""
function build_masked_laplacian(nx::Int, ny::Int, dx::Float64, dy::Float64, cyto_mask::Vector{Bool})
    # Initialize sparse matrix
    N = nx * ny
    rows = Int[]
    cols = Int[]
    vals = Float64[]
    
    # Reshape mask to 2D for easier neighbor checking
    mask_2d = reshape(cyto_mask, ny, nx)
    
    # Build Laplacian matrix
    for j in 1:ny
        for i in 1:nx
            idx = (j-1)*nx + i
            
            # Skip if point is not in cytoplasm
            if !mask_2d[j,i]
                # For points outside cytoplasm, set diagonal to 1.0 to maintain matrix structure
                push!(rows, idx)
                push!(cols, idx)
                push!(vals, 1.0)
                continue
            end
            
            # For points in cytoplasm, compute diffusion terms
            diag_val = 0.0
            
            # Left neighbor
            if i > 1
                if mask_2d[j,i-1]  # If neighbor is also in cytoplasm
                    push!(rows, idx)
                    push!(cols, idx-1)
                    push!(vals, 1.0/dx^2)
                    diag_val -= 1.0/dx^2
                end
            end
            
            # Right neighbor
            if i < nx
                if mask_2d[j,i+1]  # If neighbor is also in cytoplasm
                    push!(rows, idx)
                    push!(cols, idx+1)
                    push!(vals, 1.0/dx^2)
                    diag_val -= 1.0/dx^2
                end
            end
            
            # Bottom neighbor
            if j > 1
                if mask_2d[j-1,i]  # If neighbor is also in cytoplasm
                    push!(rows, idx)
                    push!(cols, idx-nx)
                    push!(vals, 1.0/dy^2)
                    diag_val -= 1.0/dy^2
                end
            end
            
            # Top neighbor
            if j < ny
                if mask_2d[j+1,i]  # If neighbor is also in cytoplasm
                    push!(rows, idx)
                    push!(cols, idx+nx)
                    push!(vals, 1.0/dy^2)
                    diag_val -= 1.0/dy^2
                end
            end
            
            # Add diagonal entry for the center point
            push!(rows, idx)
            push!(cols, idx)
            push!(vals, diag_val)
        end
    end
    
    # Create sparse matrix
    L = sparse(rows, cols, vals, N, N)
    return L
end

# Create mitochondria mask
function create_circular_mask(nx, ny, center, radius; dx = 1.0, dy = 1.0)
    mask = zeros(Bool, ny, nx)
    for j in 1:ny, i in 1:nx
        x = (i-0.5)/nx  # Convert to normalized coordinates
        y = (j-0.5)/ny
        if (x - center[1])^2 + (y - center[2])^2 ≤ radius^2
            mask[j,i] = true
        end
    end
    return vec(mask)  # Convert to vector for easier use
end

# Function to calculate Hill equation for ATP-dependent gating
hill_equation(atp::Float64, Kd::Float64, n::Float64) = atp > 0 ? Kd^n / (atp^n + Kd^n) : 0.0  
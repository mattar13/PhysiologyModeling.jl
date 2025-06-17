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

# Create mitochondria mask
function create_mito_mask(nx, ny, center, radius; dx = 1.0, dy = 1.0)
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
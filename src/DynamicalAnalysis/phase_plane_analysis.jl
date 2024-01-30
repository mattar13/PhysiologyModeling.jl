function meshgrid(x_range::AbstractVector, y_range::AbstractVector)
    # Create matrices for X and Y coordinates
    X = repeat(reshape(x_range, 1, :), length(y_range), 1)
    Y = repeat(reshape(y_range, :, 1), 1, length(x_range))
    return X, Y
end

function phase_plane(prob::SciMLBase.AbstractSciMLProblem, xmap::AbstractVector, ymap::AbstractVector, 
    x_idx = 1,  y_idx = 2,
)
    uI = prob.u0
    phase_map = zeros(size(xmap,1), size(ymap,1), 2)
    du = similar(uI)
    
    for (idx_x, x) in enumerate(xmap), (idx_y, y) in enumerate(ymap) 
        uI[[x_idx,y_idx]].= [x, y]
        prob.f(du, uI, prob.p, 0.0)
        phase_map[idx_x, idx_y, :] = du[[x_idx,y_idx]]
    end
    return phase_map
end

## To be tested
function find_nullclines(prob::SciMLBase.AbstractSciMLProblem, xmap::AbstractVector, ymap::AbstractVector, 
    x_idx = 1, y_idx = 2, tolerance = 1e-3
)
    # Generate the phase map
    phase_map = phase_plane(prob, xmap, ymap, x_idx, y_idx)

    # Initialize containers for nullclines
    x_nullcline = []
    y_nullcline = []

    for (idx_x, x) in enumerate(xmap), (idx_y, y) in enumerate(ymap)
        # Extract derivatives
        dx, dy = phase_map[idx_x, idx_y, :]

        # Check for x-nullcline (dx/dt ≈ 0)
        if abs(dx) < tolerance
            push!(x_nullcline, (x, y))
        end

        # Check for y-nullcline (dy/dt ≈ 0)
        if abs(dy) < tolerance
            push!(y_nullcline, (x, y))
        end
    end

    return x_nullcline, y_nullcline
end

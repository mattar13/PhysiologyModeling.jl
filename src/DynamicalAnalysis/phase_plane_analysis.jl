function meshgrid(x_range::AbstractVector, y_range::AbstractVector)
    # Create matrices for X and Y coordinates
    X = repeat(reshape(x_range, 1, :), length(y_range), 1)
    Y = repeat(reshape(y_range, :, 1), 1, length(x_range))
    return X, Y
end

function phase_plane(prob::SciMLBase.AbstractSciMLProblem, xmap::AbstractVector, ymap::AbstractVector; 
    x_idx = 2,  y_idx = 3,
)
    phase_map = zeros(size(xmap,1), size(ymap,1), 2)
    
    #Define this function inside of a function to test the phase plane
    function test_prob(x, y)
        uI = prob.u0
        du = similar(uI)
        uI[[x_idx,y_idx]].= [x, y]
        prob.f(du, uI, prob.p, 0.0)
        return du[[x_idx, y_idx]]
   end

    for (idx_x, x) in enumerate(xmap), (idx_y, y) in enumerate(ymap) 
        phase_map[idx_x, idx_y, :] = test_prob(x, y)
    end
    return phase_map
end

## In order to do this, we probably need to have xmap and ymap
function find_nullclines(prob::SciMLBase.AbstractSciMLProblem, xmap::AbstractVector, ymap::AbstractVector; 
    x_idx = 2, y_idx = 3, xic = -90.0, yic = 0.0
)

    # Initialize containers for nullclines
    ysol = fill(NaN, length(xmap))
    for (ix, x) in enumerate(xmap)
        function x_test_prob(y)
            uI = prob.u0
            du = similar(uI)
            uI[x_idx] = x
            uI[[y_idx]] .= y
            prob.f(du, uI, prob.p, 0.0)
            return du[x_idx]
        end
        #Solve the test problem
        nlsol = nlsolve(x_test_prob, [yic])
        if nlsol.f_converged
            ysol[ix] = nlsol.zero[1]
        end
    end
    
    xsol = fill(NaN, length(ymap))
    for (iy, y) in enumerate(ymap)
        function y_test_prob(x)
            uI = prob.u0
            du = similar(uI)
            uI[[x_idx]] .= x
            uI[y_idx] = y
            prob.f(du, uI, prob.p, 0.0)
            return du[y_idx]
        end
        #Solve the test problem
        nlsol = nlsolve(y_test_prob, [xic])
        if nlsol.f_converged
                xsol[iy] = nlsol.zero[1]
        end
    end

    return xsol, ysol
end
using Revise
using ElectroPhysiology, PhysiologyModeling
using Pkg; Pkg.activate("test") #Activate the testing environment
using PhysiologyPlotting, GLMakie
#This uses all of the same auxillary equations

#Make 3 points

cells = zeros(3, 2)


using GLMakie

function create_branches(origin, angle, radius, depth, branches=2)
    if depth == 0
        return []
    end
    
    branch_lines = []
    for i in 1:branches
        # Calculate the angle for each branch
        branch_angle = angle + π / 2^depth * (-1)^(i-1)
        # Calculate the end point of the branch
        end_x = origin[1] + cos(branch_angle) * radius
        end_y = origin[2] + sin(branch_angle) * radius
        # Add the branch line to the list
        push!(branch_lines, ((origin[1], end_x), (origin[2], end_y)))
        # Recursively generate the branches of the next layer
        append!(branch_lines, create_branches((end_x, end_y), branch_angle, radius / 2, depth - 1, branches))
    end
    
    return branch_lines
end

function plot_radial_tree(branch_lines)
    fig = Figure()
    ax = Axis(fig[1, 1])
    
    for line in branch_lines
        lines!(ax, line[1], line[2], color=:black)
        scatter!(ax, [line[1][1], line[1][2]], [line[2][1], line[2][2]], color=:red)
    end
    
    axislegend(ax)
    ax.aspect = DataAspect()
    return fig
end

# Example usage
layers = 4
inner_spokes = 4
radius = 1.0
branches = 2
spoke_angles = LinRange(0, 2 * π, inner_spokes)
branch_lines = []

for angle in spoke_angles
    append!(branch_lines, create_branches((0.0, 0.0), angle, radius, layers, branches))
end

fig = plot_radial_tree(branch_lines)
display(fig)

import Base.size

mutable struct CellMap{T}
	xs::Vector{T}
	ys::Vector{T}
	connections::SparseMatrixCSC{T, Int64}
     strength::SparseMatrixCSC{T, Int64}
     strength_out::Vector{T}
end

# [Map generating functions] ________________________________________________________________________________________________________________________________________#
"""
     create_even_map([kwargs...])

This function creates a evenly spaced map of cells
"""
function create_even_map(;
     xmin=0.0, dx=0.1, xmax=1.0, 
     ymin=0.0, dy=0.1, ymax=1.0)
     cells = vec([(i, j) for i in xmin:dx:xmax, j in ymin:dy:ymax])
     xs = map(c->c[1], cells)
     ys = map(c->c[2], cells)
     xs, ys
end

function create_random_map(n_cells;      
     xmin = 0.0, dx = 0.05, xmax = 2.0, 
     ymin = 0.0, dy = 0.05, ymax = 2.0
)
     xs = rand(xmin:dx:xmax, n_cells)
     ys = rand(ymin:dy:ymax, n_cells)
     return xs, ys
end

function create_ring_map(n ;center = [0.0, 0.0], r = 0.05)
     cx, cy = center
     coordinates = zeros(n, 2)
     
     for i in 1:n
         angle = 2 * π * (i / n)
         x = cx + r * cos(angle)
         y = cy + r * sin(angle)
         coordinates[i, :] .= [x, y]
     end
     
     return coordinates
end

function calculate_dendrogram_distance(xs::Vector{T}, ys::Vector{T}, old_connection_list) where T <: Real
     connection_list = Tuple[]
     for connection in old_connection_list
          x_idx, y_idx = connection
          xi = xs[x_idx]
          yi = ys[x_idx]
          xj = xs[y_idx]
          yj = ys[y_idx]
          dist = euclidean_distance([xi, yi], [xj, yj])
          push!(connection_list, (x_idx, y_idx, dist))
     end
     connection_list
end

function create_dendrogram_map(radial_lines, branches, layers; 
     origin = (0.0, 0.0), radius = 0.05/layers, branch_distance = 0.65, 
     symmetric = true #this should be true if the connections can feedback
)
     #determine angles for the inner spokes
     x0, y0 = origin
     branch_xs = Float64[origin[1],]
     branch_ys = Float64[origin[2],]
     connections = Tuple[] 
     point_index = 2
     #Next we need to increase for layers
     parent_indices = [1]
     children_indices = []
     for layer in 1:layers
          #println("Layer $layer")
          rand_branches = 0#rand(0:4)
          parent_connections = repeat(parent_indices, inner = branches+rand_branches)
          for radial_line in 1:radial_lines
               angle = 2 * π * (radial_line / radial_lines) # Start from the root
               if layer == 1
                    x = x0 + radius * cos(angle) * layer
                    y = y0 + radius * sin(angle) * layer
                    push!(branch_xs, round(x, digits = 5))
                    push!(branch_ys, round(y, digits = 5))
                    push!(connections, (parent_indices[1], point_index)) #connect all the 
                    if symmetric
                         push!(connections, (point_index, parent_indices[1]))
                    end
                    push!(children_indices, point_index)
                    point_index += 1
               else
                    n_children = branches^(layer-1) #add a random factor 
                    n_children += rand_branches
                    #println("Number of branches: $n_children")
                    #println("Number of children $(length(parent_connections))")
                    for branch in LinRange(-branch_distance, branch_distance, n_children)
                         #println("All the parents $parent_connections")
                         x = (x0 + radius * cos(angle+branch) * layer)
                         y = (y0 + radius * sin(angle+branch) * layer)
                         push!(branch_xs, round(x, digits = 5))
                         push!(branch_ys, round(y, digits = 5))
                         parent = popfirst!(parent_connections)
                         push!(connections, (parent, point_index))
                         if symmetric
                              push!(connections, (point_index, parent))
                         end
                         push!(children_indices, point_index)
                         point_index += 1
                    end
               end
          end
          #println("parent: $parent_indices")
          #println("children: $children_indices")
          parent_indices = children_indices
          children_indices = []
     end
     connections = calculate_dendrogram_distance(branch_xs, branch_ys, connections)
     branch_xs, branch_ys, connections
end

# [Connection generation functions] ___________________________________________________________________________________________________________________________-#
function connect_neighbors_radius(xs::Vector{T}, ys::Vector{T}, radii::T; constant = nothing) where T <: Real
     connections = Tuple[]  
     n_xpoints = size(xs, 1)
     n_ypoints = size(ys, 1)
     for i in 1:n_xpoints
          cell_distances = [euclidean_distance([xs[i], ys[i]], [xs[j], ys[j]]) for j in 1:n_ypoints]
          within_radius_indices = findall(d -> d <= radii, cell_distances)
          for (idx, neighbor) in enumerate(within_radius_indices)
               if isnothing(constant)
                    push!(connections, (i, neighbor, cell_distances[idx]))
               else
                    push!(connections, (i, neighbor, constant))
               end
          end
     end
     connections
end

function connect_neighbors_radius(xs::Vector{T}, ys::Vector{T}, radii::Vector{T}; constant = nothing) where T <: Real
     connections = Tuple[]  
     n_xpoints = size(xs, 1)
     n_ypoints = size(ys, 1)
     for i in 1:n_xpoints
          cell_distances = [euclidean_distance([xs[i], ys[i]], [xs[j], ys[j]]) for j in 1:n_ypoints]
          within_radius_indices = findall(d -> d <= radii[i], cell_distances)
          #println(within_radius_indices)
          for (idx, neighbor) in enumerate(within_radius_indices)
               if isnothing(constant)
                    push!(connections, (i, neighbor, cell_distances[idx]))
               else
                    push!(connections, (i, neighbor, constant))
               end
          end
     end
     connections
end

function connection_matrix(connections_list::AbstractArray{Tuple}; m = nothing, n = nothing)
     rows = map(c -> c[2], connections_list)
     cols = map(c -> c[1], connections_list)
     data = map(c -> c[3], connections_list)
     if isnothing(m)
          m = maximum(rows)
     end
     if isnothing(n)
          n = maximum(cols)
     end
     connections = sparse(rows, cols, data, m, n)
     return connections
end

connection_matrix(xs, ys, connections) = (xs, ys, connection_matrix(connections))


# [Distance calculations] __________________________________________________________________________________________________________________#

# Calculate Euclidean distance between two points
function euclidean_distance(p1::Tuple{T,T}, p2::Tuple{T,T}) where T <: Real
     x1, y1 = p1
     x2, y2 = p2
     return sqrt((x2 - x1)^2 + (y2 - y1)^2)
end

function euclidean_distance(p1, p2)
     return sqrt(sum((p1 .- p2) .^ 2))
end

function return_connected_indices(cells, connections)
     indices = Vector{Vector{Int64}}()
     n = size(cells, 1);
     for cellx in 1:n
          connected_indices = findall(x -> x != 0, connections[:, cellx])
          push!(indices, connected_indices)

     end
     return indices
end

linDist(x; m = 1.0, b = 0.0) = m*x + b 
#This is our non-linear distance function
ring(d; density = 0.05, max_dist = 0.15, slope = 0.01) = density * exp(-((d - max_dist)^2) / (2 * slope^2))

function circle_overlap_area(d, r1, r2)
     if d >= r1 + r2
         return 0.0  # No overlap
     elseif d == 0.0
          return 0.0 #Cells are one in the same
     elseif d <= abs(r1 - r2)
          return π * min(r1, r2)^2  # One circle is completely inside the other
     else
         # Calculate overlap area
         part1 = r1^2 * acos((d^2 + r1^2 - r2^2) / (2 * d * r1))
         part2 = r2^2 * acos((d^2 + r2^2 - r1^2) / (2 * d * r2))
         part3 = 0.5 * sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
         return part1 + part2 - part3
     end
end

function ring_circle_overlap_area(d; density = 1.0, r_inner = 0.1, r_outer = 0.2, r_circle = 0.2 )
     # Area of overlap between outer circle of the ring and the circle
     outer_overlap = circle_overlap_area(d, r_outer, r_circle)

     # Area of overlap between inner circle of the ring and the circle
     inner_overlap = circle_overlap_area(d, r_inner, r_circle)

     # The overlapping area with the ring is the difference
     return density*max(0.0, outer_overlap - inner_overlap)
end

function ring_circle_overlap_area(p1::Tuple{T, T}, p2::Tuple{T, T}; density = 1.0, r_inner = 0.05, r_outer = 0.09, r_circle = 0.09) where T <: Real
     # Area of overlap between outer circle of the ring and the circle
     d = euclidean_distance(p1, p2)
     outer_overlap = circle_overlap_area(d, r_outer, r_circle)

     # Area of overlap between inner circle of the ring and the circle
     inner_overlap = circle_overlap_area(d, r_inner, r_circle)

     # The overlapping area with the ring is the difference
     return density*max(0.0, outer_overlap - inner_overlap)
end

# [Constructor functions] _____________________________________________________________________________________________________________________________#
function CellMap(xs::Vector{T}, ys::Vector{T}, connections::SparseMatrixCSC{T, Int64}; 
     distance_function = ring_circle_overlap_area
) where T <: Real

     #Determine the strength of the connection via a distance function
     rows, cols, values = findnz(connections)
     new_values = map(x -> distance_function(x), values)
     strength = sparse(rows, cols, new_values, length(xs), length(ys))
     strength_out = -sum(strength, dims=1) |> vec #should we do dims 1 or dims 2

     return CellMap(xs, ys, connections, strength, strength_out)
end

# [Some utility functions for CellMap] ____________________________________________________________________________________________________________#
size(cell_map::CellMap) = length(cell_map.xs) 

function map_points(cell_map::CellMap)
     hcat(cell_map.xs, cell_map.ys)
end

function rasterize(map::CellMap; dx = 0.2, dy = 0.2)
     #round each xs and ys to integers. Then color the heatmap in based on that
     domains = map.domains
     x_map = domains[:x][1]:dx:domains[:x][2]
     y_map = domains[:y][1]:dy:domains[:y][2]
     nx = length(x_map)
     ny = length(y_map)
     grid = zeros(nx, ny)
     for i in 1:length(map.xs)
          x = round(Int64, ((map.xs[i] - domains[:x][1])/(domains[:x][2] - domains[:x][1]))/dx+1)
          y = round(Int64, ((map.ys[i] - domains[:y][1])/(domains[:y][2] - domains[:y][1]))/dy+1)
          grid[x, y] += map.vs[i]
     end
     
     return x_map, y_map, grid
end
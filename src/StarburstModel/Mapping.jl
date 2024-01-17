import Base.size
# Generate random points
function even_map(;xmin=0.0, dx=0.1, xmax=1.0, ymin=0.0, dy=0.1, ymax=1.0)
     cells = vec([(i, j) for i in xmin:dx:xmax, j in ymin:dy:ymax])
     xs = map(c->c[1], cells)
     ys = map(c->c[2], cells)
     return  [xs ys]
end

function random_map(n_points)
     return rand(n_points, 2) 
end

# Calculate Euclidean distance between two points
function euclidean_distance(p1::Tuple{T,T}, p2::Tuple{T,T}) where T <: Real
     x1, y1 = p1
     x2, y2 = p2
     return sqrt((x2 - x1)^2 + (y2 - y1)^2)
end

function euclidean_distance(p1, p2)
     return sqrt(sum((p1 .- p2) .^ 2))
end

"""
Will return either 1 for if the point has one edge, or 2 if it has 2
"""
function isboundary(point::Tuple{T, T}, domain_x::Tuple{T, T}, domain_y::Tuple{T, T}) where T <: Real
	x, y = point
	bx = Int64(domain_x[1] == x || x == domain_x[2])
	by = Int64(domain_y[1] == y || y == domain_y[2])
	return bx + by
end

function find_neighbors_radius(points, radii)
     n_points = size(points, 1)
     neighbors = Vector{Vector{Int}}(undef, n_points)
     for i in 1:n_points
          cell_distances = [euclidean_distance(points[i, :], points[j, :]) for j in 1:n_points]
          within_radius_indices = findall(d -> d <= radii[i], cell_distances)
          neighbors[i] = within_radius_indices
     end
     neighbors
end

# Create a sparse matrix to represent the connections
function create_sparse_matrix(neighbors::Vector{Vector{Int64}}, points::Matrix{T}) where T <: Real
     n_points = length(neighbors)
     I = Int[]
     J = Int[]
     V = Float64[]
     for i in 1:n_points
         for j in neighbors[i]
             push!(I, i)
             push!(J, j)
             push!(V, euclidean_distance(points[i, :], points[j, :]))
         end
     end
     return sparse(I, J, V)
end

function create_sparse_matrix(neighbors::Vector{Vector{Int64}}, points::Vector{Tuple{T,T}}) where T <: Real
     n_points = length(neighbors)
     I = Int[]
     J = Int[]
     V = Float64[]
     for i in 1:n_points
          for j in neighbors[i]
               push!(I, i)
               push!(J, j)
               push!(V, euclidean_distance(points[i], points[j]))
          end
     end
     return sparse(I, J, V)
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

mutable struct CellMap{T}
	xs::Vector{T}
	ys::Vector{T}
     radius::Vector{T}
	connections::SparseMatrixCSC{T, Int64}
     strength::SparseMatrixCSC{T, Int64}
     strength_out::Vector{T}
     domain_x::Tuple{T, T}
     domain_y::Tuple{T, T}
end

function Make_GPU!(MAP::CellMap{T}) where T <: Real



end
#These are the distance functions we can use to calculate the strength
#This is our non-linear distance function
ring(d; max_strength = 0.05, max_dist = 0.15, slope = 0.01) = max_strength * exp(-((d - max_dist)^2) / (2 * slope^2))


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

function ring_circle_overlap_area(d; 
     density = 1.0, 
     r_inner = 0.085, r_outer = 0.09, r_circle = 0.09
)
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

function CellMap(cells::Matrix{T}, radii::Vector{T}; 
     distance_function = ring_circle_overlap_area,
     domain_x = (0.0, 1.0), domain_y = (0.0, 1.0)
) where T <: Real
     neighbors = find_neighbors_radius(cells, radii)    
     connections = create_sparse_matrix(neighbors, cells)
     dropzeros!(connections)

     #Determine the strength of the connection via a distance function
     rows, cols, values = findnz(connections)

     new_values = map(x -> distance_function(x), values)
     strength = sparse(rows, cols, new_values)
     strength_out = -sum(strength, dims=2) |> vec #Preallocate the diffusion out for more efficient calculations

     return CellMap(cells[:, 1], cells[:,2], radii, connections, strength, strength_out, domain_x, domain_y)
end

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

#Some utility functions for CellMap
size(cell_map::CellMap) = length(cell_map.xs)
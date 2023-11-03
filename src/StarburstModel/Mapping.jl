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

mutable struct CellMap{T}
	xs::Vector{T}
	ys::Vector{T}
     radius::Vector{T}
	connections::SparseMatrixCSC{T, Int64}
     strength::SparseMatrixCSC{T, Int64}
     domain_x::Tuple{T, T}
     domain_y::Tuple{T, T}
end

#This is our non-linear distance function
δX(x, a, b, c) = a * exp(-((x - b)^2) / (2 * c^2))

function CellMap(cells::Matrix{T}, radii::Vector{T}; 
     max_strength = 0.05, max_dist = 0.15, slope_strength = 0.01, 
     domain_x = (0.0, 1.0), domain_y = (0.0, 1.0)
) where T <: Real
     neighbors = find_neighbors_radius(cells, radii)    
     connections = create_sparse_matrix(neighbors, cells)
     dropzeros!(connections)
     #Determine the strength of the connection via a distance function
     rows, cols, values = findnz(connections)
     println(length(rows))
     new_values = map(x -> δX(x, max_strength, max_dist, slope_strength), values)
     strengths = sparse(rows, cols, new_values)

     return CellMap(cells[:, 1], cells[:,2], radii, connections, strengths, domain_x, domain_y)
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
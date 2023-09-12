# Generate random points
function generate_points(n_points; method = :randrng,  
     xmin = -1.0, dx = 0.00001, xmax = 1.0, 
     ymin = -1.0, dy = 0.00001, ymax = 1.0
)

     if method == :rand
          points = rand(n_points, 2)     
     elseif method == :randrng
          points = [rand(xmin:dx:xmax, n_points) rand(ymin:dy:ymax, n_points)]
     elseif method == :even
          points = zeros(n_points, 2)
          n_side = Int(sqrt(n_points))
          xs = LinRange(xmin, xmax, n_side)
          ys = LinRange(ymin, ymax, n_side)
          idx = 0
          for x in xs, y in ys
               idx += 1
               points[idx, 1] = x
               points[idx, 2] = y
          end
     elseif method == :jitter
          jitter = [rand(xmin:dx:xmax, n_points) rand(ymin:dy:ymax, n_points)]
          points = zeros(n_points, 2)
          n_side = Int(sqrt(n_points))
          xs = LinRange(xmin, xmax, n_side)
          ys = LinRange(ymin, ymax, n_side)
          idx = 0
          for x in xs, y in ys
               idx += 1
               points[idx, 1] = x
               points[idx, 2] = y
          end
          points = points .+ jitter
     end
     #add boundary points
     #points = vcat(points, [xmin ymin])
     #points = vcat(points, [xmax ymax])
     #points = vcat(points, [xmin ymax])
     #points = vcat(points, [xmax ymin])
     return points
end

 # Calculate Euclidean distance between two points
function euclidean_distance(p1, p2)
     return sqrt(sum((p1 .- p2) .^ 2))
end

# Find the 5 nearest neighbors for each point
function find_nearest_neighbors(points, k)
     n_points = size(points, 1)
     neighbors = Vector{Vector{Int}}(undef, n_points)
     for i in 1:n_points
         distances = [euclidean_distance(points[i, :], points[j, :]) for j in 1:n_points]
         sorted_indices = sortperm(distances)
         neighbors[i] = sorted_indices[2:k+1]
     end
     return neighbors
end

 # Create a sparse matrix to represent the connections
function create_sparse_matrix(neighbors, points)
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

mutable struct CellMap{T}
	xs::Vector{T}
	ys::Vector{T}
	connections::SparseMatrixCSC{T, Int64}
     strength::SparseMatrixCSC{T, Int64}
     domains::NamedTuple{(:x, :y, :dx, :dy), Tuple{Tuple{T, T}, Tuple{T, T}, T, T}}
end

function CellMap(n_cells, k; 
     rand_strength_scale = 1.0, strength = 0.05, 
     map_spacing = :even,
     x = (-1.0, 1.0), y = (-1.0, 1.0), dx = 0.1, dy = 0.1
)
     domains = (x = x, y = y, dx = dx, dy = dy)     
     cells = generate_points(n_cells; 
          xmin = domains[:x][1], xmax = domains[:x][2], dx = domains[:dx],
          ymin = domains[:y][1], ymax = domains[:y][2], dy = domains[:dy],
          method = map_spacing)
     neighbors = find_nearest_neighbors(cells, k)
     
     connections = create_sparse_matrix(neighbors, cells)
     mx, my = size(connections)
     if strength == :rand
          strengths = (connections .> 0.0) .* (1 .- rand(mx,my)*rand_strength_scale)
     elseif isa(strength, Real)
          strengths = (connections .> 0.0) .* strength
     elseif strength == :ones
          strengths = (connections .> 0.0) .* 1.0
     end

     return CellMap(cells[:, 1], cells[:,2], connections, strengths, domains)
end

function ∇α(du, u, cell_map::CellMap)
     for cellx in eachindex(u)
          connected_indices = findall(x -> x != 0, cell_map.connections[:, cellx])
          flow_out = -(cell_map.strength[:, cellx].*cell_map.connections[:, cellx]) * u[cellx]
          du[cellx] += sum(flow_out) #+ sum(recieve_u)
          for (i, celly) in enumerate(connected_indices)
               flow_in = cell_map.strength[celly, cellx]*cell_map.connections[celly, cellx] * u[cellx]
               du[celly] += flow_in
          end
     end
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

function SAC_eqs_PDE(du, u, p, t, fODE, Cell_Map)
	for x in axes(u, 2)
		du_i = du[:, x]
		u_i = u[:, x]
		fODE(du_i, u_i, p, t)
		du[:, x] .= du_i 
	end
     #Conduct diffusion for the correct terme
     dE = du[6, :]
     println(dE)
     E = u[6, :]
     ∇α(dE, E, Cell_Map)
     du[6, :] .= dE
end
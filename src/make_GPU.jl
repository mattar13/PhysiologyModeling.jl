function make_GPU(MAP::CellMap{T}) where T <: Real
	xs = MAP.xs |> Vector{Float32}
     ys = MAP.ys |> Vector{Float32}
     radius = MAP.radius |> Vector{Float32}
     connections = MAP.connections .|> Float32 |> CUSPARSE.CuSparseMatrixCSC
     strength = MAP.strength .|> Float32 |> CUSPARSE.CuSparseMatrixCSC
     strength_out = MAP.strength_out |> Vector{Float32}
     return CellMap(xs, ys, radius, connections, strength, strength_out)
end
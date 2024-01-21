mutable struct CellMapGPU{T}
	xs::Vector{T}
	ys::Vector{T}
     radius::Vector{T}
	connections::CUDA.CUSPARSE.CuSparseMatrixCSC{T, Int32}
     strength::CUDA.CUSPARSE.CuSparseMatrixCSC{T, Int32}
     strength_out::Vector{T}
end

function make_GPU(MAP::CellMap{T}) where T <: Real
	xs = MAP.xs |> Vector{Float32}
     ys = MAP.ys |> Vector{Float32}
     radius = MAP.radius |> Vector{Float32}
     connections = MAP.connections .|> Float32 |> CUSPARSE.CuSparseMatrixCSC
     strength = MAP.strength .|> Float32 |> CUSPARSE.CuSparseMatrixCSC
     strength_out = MAP.strength_out |> Vector{Float32}
     return CellMapGPU(xs, ys, radius, connections, strength, strength_out)
end

function SAC_PDE_GPU(du, u, p, t, MAP::CellMapGPU)
     CUDA.@sync SAC_PDE(du, u, p, t, MAP)
     return
end
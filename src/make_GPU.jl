mutable struct CellMapGPU{T}
	xs::Vector{T}
	ys::Vector{T}
	connections::CUDA.CUSPARSE.CuSparseMatrixCSC{T, Int32}
     strength::CUDA.CUSPARSE.CuSparseMatrixCSC{T, Int32}
     strength_out::CuArray{T}
end

function make_GPU(MAP::CellMap{T}) where T <: Real
	xs = MAP.xs |> Vector{Float32}
     ys = MAP.ys |> Vector{Float32}
     connections = MAP.connections .|> Float32 |> CUSPARSE.CuSparseMatrixCSC
     strength = MAP.strength .|> Float32 |> CUSPARSE.CuSparseMatrixCSC
     strength_out = MAP.strength_out |> CuArray{Float32}
     return CellMapGPU(xs, ys, connections, strength, strength_out)
end

function DIFFUSION_MODEL_GPU(du, u, p, t; active_cell = 221)
     CUDA.@sync DIFFUSION_MODEL(du, u, p, t; active_cell = active_cell)
     #CUDA.@sync ∇α(du, u, p, t)
     return
end

function SAC_PDE_GPU(du, u, p, t, MAP::CellMapGPU)
     CUDA.@sync SAC_PDE(du, u, p, t, MAP)
     return
end
using BenchmarkTools
using CUDA
using Test

N = 2^20
x = fill(1.0f0, N)  # a vector filled with 1.0 (Float32)
y = fill(2.0f0, N)  # a vector filled with 2.0

x_d = CUDA.fill(1.0f0, N)  # a vector stored on the GPU filled with 1.0 (Float32)
y_d = CUDA.fill(2.0f0, N)  # a vector stored on the GPU filled with 2.0

function sequential_add!(y, x)
    for i in eachindex(y, x)
        @inbounds y[i] += x[i]
    end    
    return nothing
end    

function parallel_add!(y, x)
    Threads.@threads for i in eachindex(y, x)
        @inbounds y[i] += x[i]
    end    
    return nothing
end    

function add_broadcast!(y, x) #Using this method seems to be the fastest way to do operations
    CUDA.@sync y .+= x
    return
end

function gpu_add1!(y, x)
    for i in eachindex(y)
        @inbounds y[i] += x[i]
    end
    return nothing
end

function bench_gpu1!(y, x)
    CUDA.@sync begin
        @cuda gpu_add1!(y, x)
    end
end

@cuda gpu_add1!(y_d, x_d)

@btime sequential_add!($y, $x)
@btime parallel_add!($y, $x)
@btime add_broadcast!($y_d, $x_d)
@btime bench_gpu1!($y_d, $x_d)
nothing
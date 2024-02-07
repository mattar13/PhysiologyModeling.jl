
#This is a useful function for expanding the 
function extend_dims(A,which_dim)
     s = [size(A)...]
     insert!(s,which_dim,1)
     return reshape(A, s...)
end

import ElectroPhysiology.Experiment
function Experiment(sol::SciMLBase.AbstractSciMLSolution; 
     channels = :vars,
)
     sol_arr = extend_dims(sol |> Array, 3)
     if channels == :vars
          sol_arr = permutedims(sol_arr, (3,2,1))
     end
     #Enter the Header here
     Experiment(sol.t, sol_arr)
end

extract_dict(d::Dict{String, Float64}, keys) = map(k -> d[k], keys)
extract_p0(d::Dict{String, Float64})  = extract_dict(d, keys_p0)
extract_p0(d::NamedTuple) = d
extract_p0(d) = d

extract_u0(d::Dict{String, Float64}) = extract_dict(d, keys_u0)

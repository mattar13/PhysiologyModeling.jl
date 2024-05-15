#This is a useful function for expanding the 
function extend_dims(A,which_dim)
     s = [size(A)...]
     insert!(s,which_dim,1)
     return reshape(A, s...)
end

import ElectroPhysiology.Experiment
function Experiment(sol::SciMLBase.AbstractSciMLSolution;
     dt = nothing, 
     channels = :vars,
)    
     if !isnothing(dt)
          model_t = sol.t[1]:dt:sol.t[end] |> collect
          new_sol = sol.(model_t)
          model_arr = reshape(hcat(new_sol...), 1, size(model_t,1), size(sol, 1))
     else
          model_t = sol.t
          model_arr = sol
          model_arr = extend_dims(model_arr |> Array, 3)
     end
     if channels == :vars
          model_arr = permutedims(model_arr, (3,2,1))
     end
     #Enter the Header here
     Experiment(model_t, model_arr)
end

#This is a useful function for expanding the 
function extend_dims(A,which_dim)
     s = [size(A)...]
     insert!(s,which_dim,1)
     return reshape(A, s...)
end

import ElectroPhysiology.Experiment
function Experiment(sol::SciMLBase.AbstractSciMLSolution)
     sol_arr = extend_dims(sol |> Array, 3)
     #Enter the Header here
     Experiment(sol.t, sol_arr)
end


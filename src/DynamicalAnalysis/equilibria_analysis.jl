struct equilibria_object{T}
    stable::Array{Array{T}}
    unstable::Array{Array{T}}
    saddle::Array{Array{T}}
    stable_focus::Array{Array{T}}
    unstable_focus::Array{Array{T}}
end

import Base.getindex
function getindex(eq::equilibria_object, sym::Symbol)
    if !isnothing(sym |> u_find)
        var = sym |> u_find
        a = [
            !isempty(eq.stable) ? eq.stable[1][var] : NaN,
            !isempty(eq.unstable) ? eq.unstable[1][var] : NaN,
            !isempty(eq.saddle) ? eq.saddle[1][var] : NaN,
            !isempty(eq.stable_focus) ? eq.stable_focus[1][var] : NaN,
            !isempty(eq.unstable_focus) ? eq.unstable_focus[1][var] : NaN
        ]
        return a
    elseif sym == :stable
        return eq.stable
    elseif sym == :unstable
        return eq.unstable
    elseif sym == :saddle
        return eq.saddle
    elseif sym == :stable_focus
        return eq.stable_focus
    elseif sym == :unstable_focus
        return eq.unstable_focus
    end
end

getindex(eq::equilibria_object, syms...) = map(sym -> eq[sym], syms)

function getindex(eq::equilibria_object, idx::Int64)  
    if idx == 1
        return eq.stable
    elseif idx == 2
        return eq.unstable
    elseif idx == 3
        return eq.saddle
    elseif idx == 4
        return eq.stable_focus
    elseif idx == 5
        return eq.unstable_focus
    end
end

import Base.length, Base.print
#We can import a function that counts all of the equilibria events in the object
length(eq::equilibria_object) = length(eq.stable) + length(eq.unstable) + length(eq.saddle) + length(eq.unstable_focus) + length(eq.stable_focus)

#This function displays information about the equilibrium
function print(eq::equilibria_object; msg = "Existant Equilibrium", vars = [:v])
    println(msg)
    #var_idxs = map(vr -> findall(x -> x==vr, tar_conds)[1], vars)
    var_idxs = map(vr -> vr |> u_find, vars)
    eq.unstable != [] ? println("Unstable Equilibrium: $(eq.unstable[1][var_idxs])") : nothing
    eq.stable != [] ? println("Stable Equilibrium: $(eq.stable[1][var_idxs])") : nothing
    eq.saddle != [] ? println("Saddle Equilibrium: $(eq.saddle[1][var_idxs])") : nothing
    eq.unstable_focus != [] ? println("Unstable Focus Equilibrium: $(eq.unstable_focus[1][var_idxs])") : nothing
    eq.stable_focus != [] ? println("Stable Focus Equilibrium: $(eq.stable_focus[1][var_idxs])") : nothing
end

export print, length

#Conduct a stability analysis of the current
#We have to ensure the points we pick are realistic
function find_fixed_points(prob::SciMLBase.AbstractSciMLProblem, xmap::AbstractVector, ymap::AbstractVector; 
    x_idx = 2, y_idx = 3, 
    precision = 3,
    verbose = false
    )
    xlims = (xmap[1], xmap[end])
    ylims = (ymap[1], ymap[end])
    fixed_points = Vector[] #This stores all of the fixed points
    function test_prob(ux)
        du = similar(ux)
        prob.f(du, ux, prob.p, 0.0) #solve this problem for zero
        return du       
    end
    for (idx_x, x) in enumerate(xmap), (idx_y, y) in enumerate(ymap)
        #Iterate through the y range looking for stable points
        uI = copy(prob.u0)
        uI[[x_idx, y_idx]] .= [x, y]
        if verbose
            println("Checking parameter")
            println("var 1 = $x")
            println("var 2 = $y")
        end
        res = nlsolve(test_prob, uI)
        #Check to see if the nlsolve failed
        
        #If any of the results are out of bounds then we have to cut it off
        if verbose
            print("Results of the nlsolve for uI: ")
            #println(res)
            print("f = 0: ")
            println(res.zero)
        end
        var1_inbounds = xlims[1] < res.zero[x_idx] <= xlims[2]
        var2_inbounds = ylims[1] < res.zero[y_idx] <= ylims[2]
        if var1_inbounds && var2_inbounds
            fixed_points_check = map(x -> round(x, digits=precision), res.zero)
            if fixed_points_check ∉ fixed_points && res.f_converged
                if verbose
                    println(res.f_converged)
                    println(fixed_points_check)
                end
                push!(fixed_points, fixed_points_check)
            end
                
        else
            if !var1_inbounds && verbose
                println("Variable $(x_idx) is out of bounds at $(res.zero[1]) will be skipped")
            end

            if !var2_inbounds && verbose
                println("Variable $(y_idx) is out of bounds at $(res.zero[2]) will be skipped")
            end
        end
    end
    return fixed_points
end

function find_equilibria(prob::SciMLBase.AbstractSciMLProblem, fixed_point; x_idx = 2, y_idx = 3)
    function test_prob(ux)
        du = similar(ux)
        prob.f(du, ux, prob.p, 0.0) #solve this problem for zero
        return du       
    end
    test_jac = ForwardDiff.jacobian(test_prob, fixed_point)[[x_idx, y_idx], [x_idx, y_idx]]
    ev = eigvals(test_jac)
    if sign(real(ev[1])) != sign(real(ev[2])) #If the signs of the eigen values are not equal, the equilbria is unstable   
        return :Saddle
    elseif imag(ev[1]) != 0.0 #This detects whether the equilibria is a focus equilibria
        if real(ev[1]) >= 0.0
            return :UnstableFocus
        elseif real(ev[1]) < 0.0
            return :StableFocus
        end
    else #If it is not a saddle and not a focus, then it is just a plain equilibria
        if ev[1] >= 0.0
            return :Unstable
        elseif ev[1] < 0.0
            return :Stable
        end
    end 

end
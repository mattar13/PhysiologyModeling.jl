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
function find_fixed_points(prob::SciMLBase.AbstractSciMLProblem, xmap::AbstractVector, ymap::AbstractVector; x_idx = 2, y_idx = 3, verbose = true)
    xlims = (xmap[1], xmap[end])
    ylims = (ymap[1], ymap[end])
    fixed_points = Tuple
    function test_prob(x)
        uI = prob.u0
        du = similar(uI)
        uI[[x_idx,y_idx]].= x
        prob.f(du, uI, prob.p, 0.0)
        return du[[x_idx, y_idx]]
    end
    for (idx_x, x) in enumerate(xmap), (idx_y, y) in enumerate(ymap)
        #Iterate through the y range looking for stable points
        if verbose
            println("Checking parameter")
            println("var 1 = $x")
            println("var 2 = $y")
        end
        res = nlsolve(test_prob, [x,y])
        #If any of the results are out of bounds then we have to cut it off
        if verbose
            print("Results of the nlsolve for uI: ")
            #println(res)
            print("f = 0: ")
            println(res.zero)
        end
        var1_inbounds = xlims[1] < res.zero[1] <= xlims[2]
        var2_inbounds = ylims[1] < res.zero[2] <= ylims[2]
        if var1_inbounds && var2_inbounds

        else
            if !var1_inbounds && verbose
                println("Variable $(x_idx) is out of bounds at $(res.zero[1]) will be skipped")
            end

            if !var2_inbounds && verbose
                println("Variable $(y_idx) is out of bounds at $(res.zero[2]) will be skipped")
            end
        end
    end

end


#%%
function find_equilibria(prob::SciMLBase.AbstractSciMLProblem, xlims::Tuple, ylims::Tuple;
    x_vars = 2, y_vars = 3, 
    equilibrium_resolution = 10,
    precision = 2, check_min = 1e-5, verbose = false
)
    stable = Array{Array{Float64}}([])
    unstable = Array{Array{Float64}}([])
    saddle = Array{Array{Float64}}([])
    unstable_focus = Array{Array{Float64}}([])
    stable_focus = Array{Array{Float64}}([])
    storage = Array{Array{Float64}}([])
    xmap = LinRange(xlims[1], xlims[2], equilibrium_resolution)
    ymap = LinRange(ylims[1], ylims[2], equilibrium_resolution)
    function df(ux) 
        #println("at least we made it here")
        du = similar(ux)
        prob.f(du, ux, prob.p, 0.0) #solve this problem for zero
        return du           
    end

    for (idx_x, x) in enumerate(xmap) for (idx_y, y) in enumerate(ymap)
            #Iterate through the y range looking for stable points
            if verbose
                println("Checking parameter")
                println("var 1 = $x")
                println("var 2 = $y")
            end
            uI = copy(prob.u0)
            uI[[x_vars, y_vars]] .= [x, y]
        
            if verbose
                println("Condition to check = $uI")
            end
        
            res = nlsolve(df, uI)
            #If any of the results are out of bounds then we have to cut it off
            if verbose
                print("Results of the nlsolve for uI: ")
                #println(res)
                print("f = 0: ")
                println(res.zero)
            end
        
            #We want to remove all cases where the result is out of bounds
            var1_inbounds = xlims[1] < res.zero[var_idx[1]] <= xlims[2]
            var2_inbounds = ylims[1] < res.zero[var_idx[2]] <= ylims[2]

            if var1_inbounds && var2_inbounds #We want to make sure the equilibria is not out of bounds
                equilibria = map(x -> round(x, digits=precision), res.zero)
                check = df(res.zero)
                if verbose
                    print("Checking the solution and it's proximity to zero: dUᵢ = ")
                    println(check)
                end
                if equilibria in storage
                    nothing
                elseif any(isnan, res.zero)
                    nothing
                elseif check[1] > check_min || check[2] > check_min
                    nothing
                else
                    push!(storage, equilibria)
                    #push!(eq_dict[:all], res.zero)
                    #If the equilibria is not in the stable or unstable points
                    #Test it's stability
                    J_mat = ForwardDiff.jacobian(df, res.zero)[[var_idx...], [var_idx...]]
                    ev = eigvals(J_mat)
            
                    if sign(real(ev[1])) != sign(real(ev[2]))
                        #If the sign of the eigenvalues are opposite, then it is a unstable saddle
                        #println("Saddle")
                        push!(saddle, res.zero)
                    else
                        if imag(ev[1]) != 0 #This detects whether or not the equilibria is a focus
                            if real(ev[1]) >= 0.0
                                #println("Unstable focus")
                                push!(unstable_focus, res.zero)
                            else
                                #println("Unstable focus")
                                push!(stable_focus, res.zero)
                            end
                        else
                            if ev[1] > 0.0
                                push!(unstable, res.zero)
                            else
                                push!(stable, res.zero)
                            end
                        end
                    end
                end
            else
                if !var1_inbounds && verbose
                    println("Variable $(vars[1]) is out of bounds at $(res.zero[var_idx[1]]) will be skipped")
                end

                if !var2_inbounds && verbose
                    println("Variable $(vars[2]) is out of bounds at $(res.zero[var_idx[2]]) will be skipped")
                end
            end
        end
    end
    return equilibria_object{Float64}(stable, unstable, saddle, unstable_focus, stable_focus)
end

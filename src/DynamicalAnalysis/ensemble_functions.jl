function OneVarEnsembleProb(prob, i, repeat, params)
     pI = prob.p
     pI[idx] = params[i]
     remake(prob, p = pI)
end

function InitialCond_Param_EnsembleProb(prob, i, repeat, ics, params; idx = 1)
     pI = prob.p
     pI[idx] = params[i]
     remake(prob, p = pI, u0 = ics[:, i])
end
function OneVarEnsembleProb(prob, i, repeat, parameter_range; idx = 1)
     pI = prob.p
     pI[idx] = parameter_range[i]
     remake(prob, p = pI)
end
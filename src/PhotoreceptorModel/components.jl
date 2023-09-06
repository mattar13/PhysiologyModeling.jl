@connector function Point(; name)
     sts = @variables v(t)=1.0 i(t)=1.0 [connect = Flow] #Flow sums to zero allowing continuity between branches
     ODESystem(Equation[], t, sts, []; name = name)
end

#This is a branch point of a neuron that has a voltage and current. 
#There is no resistance or conductance between the two points
@component function IdealCable(; name) 
     @named p = Point() #Positive Terminal
     @named n = Point() #Negative terminal
     sts = @variables v(t)=1.0 i(t)=1.0
     eqs = [v ~ p.v - n.v #Voltage is the voltage difference between the positive and negative terminal
          0 ~ p.i + n.i #Current on positve end has to equal curren to negative
          i ~ p.i]
     compose(ODESystem(eqs, t, sts, []; name = name), p, n)
end

#cables have resistance and capacitance
@component function Cable(; name, R = 1.0, C = 10.0)
     @named ic = IdealCable()
     @unpack v, i = ic
     ps = @parameters R = R C = C
     D = Differential(t)
     eqs = [
          D(v) ~ (R*i) / C,
     ]
     extend(ODESystem(eqs, t, [], ps; name = name), ic)
end

@component function Ground(; name) #Voltage is equal to zero on the groun
     @named g = Pin()
     eqs = [g.v ~ 0]
     compose(ODESystem(eqs, t, [], []; name = name), g)
end

@component function ConstantVoltage(; name, V = 1.0)
     @named oneport = OnePort()
     @unpack v = oneport
     ps = @parameters V = V
     eqs = [
         V ~ v,
     ]
     extend(ODESystem(eqs, t, [], ps; name = name), oneport)
end

#%% What we have built already
@component function LeakySoma(params::Dict; name, hv = 0.0, v = 0.0)
     sts = @variables hv(t)=hv v(t)=v mKV(t) = 0.430 hKV(t) = 0.999#States
     ps = @parameters cm = params[:cm] gL = params[:gL] EL = params[:EL] #Parameters
     D = Differential(t)
     #Dependant functions
     I_LEAK(v) = gL*(v-EL)
     I_KV(v) = gKV * mKV^3 * hKV * (v - EK)
     eqs = [
          D(hv) ~ -hv
          D(v) ~ -(I_LEAK(v) + hv)/cm
          D(mKV) ~ αmKV(v) * (1 - mKV) - βmKV(v) * mKV
          D(hKV) ~ αhKV(v) * (1 - hKV) - βhKV(v) * hKV 
     ]
     ODESystem(eqs, t, sts, ps; name = name)
end


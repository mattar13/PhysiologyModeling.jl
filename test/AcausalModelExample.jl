using DifferentialEquations, ModelingToolkit

#Build an acausal circuit
@variables t

#When you connect two pins, the flow variable will sum to zero
@connector function Pin(; name)
    sts = @variables v(t)=1.0 i(t)=1.0 [connect = Flow]
    ODESystem(Equation[], t, sts, []; name = name)
end

#A one port has a positive and negative terminal
@component function OnePort(; name)
     @named p = Pin() #Positive Terminal
     @named n = Pin() #Negative terminal
     sts = @variables v(t)=1.0 i(t)=1.0
     eqs = [v ~ p.v - n.v #Voltage is the voltage difference between the positive and negative terminal
            0 ~ p.i + n.i #Current on positve end has to equal curren to negative
            i ~ p.i]
     compose(ODESystem(eqs, t, sts, []; name = name), p, n)
end

#The Resistor inherits the OnePort 
@component function Resistor(; name, R = 1.0)
     @named oneport = OnePort()
     @unpack v, i = oneport #This takes the voltage and current from the oneport
     ps = @parameters R = R #An external parameter that can be set
     eqs = [
          v ~ i * R, #A resistor must satisfy Ohms law
     ]
     extend(ODESystem(eqs, t, [], ps; name = name), oneport)
end

@component function Capacitor(; name, C = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters C = C
    D = Differential(t)
    eqs = [
        D(v) ~ i / C,
    ]
    extend(ODESystem(eqs, t, [], ps; name = name), oneport)
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

R = 1.0
C = 1.0
V = 1.0
@named resistor = Resistor(R = R)
@named capacitor = Capacitor(C = C)
@named source = ConstantVoltage(V = V)
@named ground = Ground()

rc_eqs = [connect(source.p, resistor.p)
          connect(resistor.n, capacitor.p)
          connect(capacitor.n, source.n)
          connect(capacitor.n, ground.g)]

@named _rc_model = ODESystem(rc_eqs, t)
@named rc_model = compose(_rc_model,
                          [resistor, capacitor, source, ground])
sys = structural_simplify(rc_model)

u0 = [
    capacitor.v => 0.0,
]
prob = ODAEProblem(sys, u0, (0, 10.0))
prob |> typeof |> fieldnames
sys
sol = solve(prob, Tsit5())
plot(sol)
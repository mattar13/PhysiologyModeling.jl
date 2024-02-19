#= Code from https://github.com/FellerLabCodeShare/DSGC-Velocity-Project/blob/main/Fig%205%20Modeling/runModel.m
function [Vm] = runModel(gExc,gInh)
% Function to run a simple conductance model using forward euler method

% Intialize model
% Vm = NaN(paramsModel.nPts,1);
nPts = size(gInh,1);
Vm = NaN(nPts,1);
Vm(1,:) = paramsModel.eLeak;
dV = NaN;

% Integrate via forward euler method
for k = 1:(nPts - 1)%(paramsModel.nPts - 1)
    % Calculate change in V 
    dV = -( ...
        gExc(k) * (Vm(k) - paramsModel.eExc) + ... %exc conductance term
        gInh(k) * (Vm(k) - paramsModel.eInh) + ... %inh conductance term
        paramsModel.gLeak*(Vm(k) - paramsModel.eLeak)) ... % leak conductance term
        / paramsModel.cap; %capacitance
    Vm(k+1,1) = Vm(k) + dV*paramsModel.dt; %numerically integrate
end

en=#

function InhExcModel(du, u, p, t)
    (C_m, g_leak, E_leak, E_Exc, E_Inh) = p

    dv = view(du, 1)
    dgE = view(du, 2)
    dgI = view(du, 3)

    v = view(u, 1)
    gE = view(u, 2)
    gI = view(u, 3)

    @. dv = -(
        g_leak * (v - E_leak) + 
        gE * (v - E_Exc) +
        gI * (v - E_Inh)
    ) / C_m
    @. dgE = 0.0
    @. dgI = 0.0
    nothing
end


function gaussian(x; σ = 1.0, μ = 1.0)
    1/√(2*π*σ)*exp((-x-μ^2)/(2*σ^2))
end
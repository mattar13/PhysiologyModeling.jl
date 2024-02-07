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

function IncExcModel(du, u, p, t)



end

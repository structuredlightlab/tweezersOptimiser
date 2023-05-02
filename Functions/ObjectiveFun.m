function [obj,grad] = ObjectiveFun(x,Ops,Case,optmType)
%ObjectiveFun evaluates the objective function and its gradient using the GWS operators, 
%given a specific Case

% INPUT
% x             optimisation variable; it is either
%               (a) real and imaginary parts of the complex amplitudes
%               of the modes in the far field; 1D array, first half are the real and second
%               are the imaginary parts OR
%               (b) the phase of the far-field modes
% Ops           structure containing the relevant GWS operators and their transposes
% Case          string specifying objective and constraints
% optmType      what kind of optimiser do you want to run? more details in 'Optimiser_phaseOnly_FUN'

% OUTPUT
% obj           value of the objective function
% grad          gradient of the objective function

%objective function
objStr = extractBefore(Case,'_st:');

switch optmType
    case 'phase_op1'
        normF = sqrt(x(1:end/2).^2 + x(end/2+1:end).^2);
        x = x./[normF; normF];
        x = x/sqrt(length(x)/2);
    case 'phase_op3'
        x = exp(1i*x);
end
xCT = x';
xConj = conj(x);


switch objStr
    case 'obj:Kx'
       [obj, grad] = OptmForceStiffness(x,xConj,xCT,Ops.Dx,Ops.DxT,optmType);
    otherwise
        error('specified objective not supported')
end


end

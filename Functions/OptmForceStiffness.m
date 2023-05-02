function [fs,grad] = OptmForceStiffness(x,xConj,xCT,op,opT,optmType)
%OptmForceStiffness calculates the force or stiffness and their gradient
%based on the GWS operators

% INPUT
% x                  optimisation variable (processed as needed)
% xConj              conjugate of the optimisation variable
% xCT                conjugate transpose of the optimisation variable
% op                 the relevant operator (Q,K,B or D); 2D array
% opT                conjugate transpose of op

% OUTPUT
% fs                 force or stiffness (dependong on the op); scalar
% grad               gradient of fs (w.r.t. x); 1D array same size as x



%force/stiffness
fs = real(xCT*op*x); 

%gradient
switch optmType
    case {'cmplxAmpl','phase_op1','phase_op2'}
        grad = real((op + opT)*x);
    case 'phase_op3'
        grad = real(-1i*xConj.*op*x + 1i*x.*opT*xConj);
end


end

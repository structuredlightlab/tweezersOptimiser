function [cIneq,cEq,cIneqG,cEqG] = nonLinConstraints(x,Ops,Case,optmType,R,a)
%this function provides the constraint functions for the
%'OptimiserForStiffness_usingGWS_solverBased' code

% INPUT
% x             optimisation variable; it is either
%               (a) real and imaginary parts of the complex amplitudes
%               of the modes in the far field; 1D array, first half are the real and second
%               are the imaginary parts OR
%               (b) the phase of the far-field modes
% Ops           structure containing the relevant GWS operators and their transposes
% Case          string specifying objective and constraints
% optmType      what kind of optimiser do you want to run? more details in 'Optimiser_phaseOnly_FUN'
% R             desired ratio between z and x stiffness

% OUTPUT
% cIneq
% cEq           
% cIneqG        
% cEqG

%optimisation variable
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

%string specifying the constraints
conStr = extractAfter(Case,'_st:');

%---calculate the force, stiffness and gradients(of operators)
if contains(conStr,'Fx');   [fx,g_fx]     = OptmForceStiffness(x,xConj,xCT,Ops.Bx,Ops.BxT,optmType);     end
if contains(conStr,'Fy');   [fy,g_fy]     = OptmForceStiffness(x,xConj,xCT,Ops.By,Ops.ByT,optmType);     end
if contains(conStr,'Fz');   [fz,g_fz]     = OptmForceStiffness(x,xConj,xCT,Ops.Bz,Ops.BzT,optmType);     end
if contains(conStr,'Kx');   [kx,g_kx]     = OptmForceStiffness(x,xConj,xCT,Ops.Dx,Ops.DxT,optmType);     end
if contains(conStr,'Ky');   [ky,g_ky]     = OptmForceStiffness(x,xConj,xCT,Ops.Dy,Ops.DyT,optmType);     end
if contains(conStr,'Kz');   [kz,g_kz]     = OptmForceStiffness(x,xConj,xCT,Ops.Dz,Ops.DzT,optmType);     end
if contains(conStr,'Kxy');  [kxy,g_kxy]   = OptmForceStiffness(x,xConj,xCT,Ops.Dxy,Ops.DxyT,optmType);   end
if contains(conStr,'Kxny'); [kxny,g_kxny] = OptmForceStiffness(x,xConj,xCT,Ops.Dxny,Ops.DxnyT,optmType); end
if contains(conStr,'Kxz');  [kxz,g_kxz]   = OptmForceStiffness(x,xConj,xCT,Ops.Dxz,Ops.DxzT,optmType);   end
if contains(conStr,'Kxnz'); [kxnz,g_kxnz] = OptmForceStiffness(x,xConj,xCT,Ops.Dxnz,Ops.DxnzT,optmType); end
if contains(conStr,'Kyz');  [kyz,g_kyz]   = OptmForceStiffness(x,xConj,xCT,Ops.Dyz,Ops.DyzT,optmType);   end
if contains(conStr,'Kynz'); [kynz,g_kynz] = OptmForceStiffness(x,xConj,xCT,Ops.Dynz,Ops.DynzT,optmType); end

%%% INEQUALITY CONSTRAINTS %%%
cIneq = []; cIneqG = [];
if contains(conStr,'Kx<=0'); cIneq = [cIneq, kx]; cIneqG = [cIneqG, g_kx]; end


%%% EQUALITY CONSTRAINTS %%%
switch optmType
    case {'cmplxAmpl','phase_op1'}
        %keep far-field power constant
        cEq(1) = sum(x.^2) - 1;    cEqG(:,1) = 2*x;
    case 'phase_op2'
        cEq = transpose( x(1:end/2).^2 + x(end/2+1:end).^2 - a.^2 ); 
        cEqG = 2*[diag(x(1:end/2)); diag(x(end/2+1:end))];
    case 'phase_op3'
        cEq=[]; cEqG=[];
end

%other constraitns
if contains(conStr,'Ky=Kx');    cEq = [cEq, ky-kx];    cEqG = [cEqG, g_ky-g_kx];    end
if contains(conStr,'Kxy=Kx');   cEq = [cEq, kxy-kx];   cEqG = [cEqG, g_kxy-g_kx];   end
if contains(conStr,'Kxny=Kx');  cEq = [cEq, kxny-kx];  cEqG = [cEqG, g_kxny-g_kx];  end
if contains(conStr,'Kz=Kx*R');  cEq = [cEq, kz-kx*R];  cEqG = [cEqG, g_kz-g_kx*R];  end
if contains(conStr,'Kxnz=Kxz'); cEq = [cEq, kxnz-kxz]; cEqG = [cEqG, g_kxnz-g_kxz]; end
if contains(conStr,'Kynz=Kyz'); cEq = [cEq, kynz-kyz]; cEqG = [cEqG, g_kynz-g_kyz]; end
if contains(conStr,'Kynz=Kxz'); cEq = [cEq, kynz-kxz]; cEqG = [cEqG, g_kynz-g_kxz]; end
if contains(conStr,'Kyz=Kxz');  cEq = [cEq, kyz-kxz];  cEqG = [cEqG, g_kyz-g_kxz];  end
if contains(conStr,'Fx=0');     cEq = [cEq, fx];       cEqG = [cEqG, g_fx];         end
if contains(conStr,'Fy=0');     cEq = [cEq, fy];       cEqG = [cEqG, g_fy];         end
if contains(conStr,'Fz=0');     cEq = [cEq, fz];       cEqG = [cEqG, g_fz];         end


end


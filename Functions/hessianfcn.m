function [Hout] = hessianfcn(x,lambda,Ops,Case,optmType,R,Ha)

%note to self: x and lambda have to be included in the inputs even if they
%are not used in the claculation; matlab expects these two inputs

% INPUT
% x             optimisation variable; it is either
%               (a) real and imaginary parts of the complex amplitudes
%               of the modes in the far field; 1D array, first half are the real and second
%               are the imaginary parts OR
%               (b) the phase of the far-field modes
% Ops           structure containing the relevant GWS operators and their transposes
% Case          string specifying objective and constraints
% optmType      what kind of optimiser do you want to run? more details in 'Optimiser_FUN'
% R             desired ratio between z and x stiffness

% OUTPUT
% Hout          Hessian of the Lagrangian, made up from Hessian of
%               the objective function Hobj, Hessians of
%               the equality constraints Heq, and Hessians of the
%               inequality constraints Hineq; 2D array


switch optmType
    case 'phase_op1'
        normF = sqrt(x(1:end/2).^2 + x(end/2+1:end).^2);
        x = x./[normF; normF];
        x = x/sqrt(length(x)/2);
    case 'phase_op3'
        x = exp(1i*x);
end
xT = transpose(x);
xConj = conj(x);

%---calculate the Hessians(of operators)
if contains(Case,'Fx');   H_Bx   = hessianGWS(Ops.Bx,Ops.BxT,x,xConj,xT,optmType);     end
if contains(Case,'Fy');   H_By   = hessianGWS(Ops.By,Ops.ByT,x,xConj,xT,optmType);     end
if contains(Case,'Fz');   H_Bz   = hessianGWS(Ops.Bz,Ops.BzT,x,xConj,xT,optmType);     end
if contains(Case,'Kx');   H_Dx   = hessianGWS(Ops.Dx,Ops.DxT,x,xConj,xT,optmType);     end
if contains(Case,'Ky');   H_Dy   = hessianGWS(Ops.Dy,Ops.DyT,x,xConj,xT,optmType);     end
if contains(Case,'Kz');   H_Dz   = hessianGWS(Ops.Dz,Ops.DzT,x,xConj,xT,optmType);     end
if contains(Case,'Kxy');  H_Dxy  = hessianGWS(Ops.Dxy,Ops.DxyT,x,xConj,xT,optmType);   end
if contains(Case,'Kxny'); H_Dxny = hessianGWS(Ops.Dxny,Ops.DxnyT,x,xConj,xT,optmType); end
if contains(Case,'Kxz');  H_Dxz  = hessianGWS(Ops.Dxz,Ops.DxzT,x,xConj,xT,optmType);   end
if contains(Case,'Kxnz'); H_Dxnz = hessianGWS(Ops.Dxnz,Ops.DxnzT,x,xConj,xT,optmType); end
if contains(Case,'Kyz');  H_Dyz  = hessianGWS(Ops.Dyz,Ops.DyzT,x,xConj,xT,optmType);   end
if contains(Case,'Kynz'); H_Dynz = hessianGWS(Ops.Dynz,Ops.DynzT,x,xConj,xT,optmType); end

%---process the Case string to determine which Hessians (of functions) are needed
%%% OBJECTIVE %%%
objStr = extractBefore(Case,'_st:');
switch objStr
    case 'obj:Kx'
        Hobj = H_Dx;
    otherwise
        error('specified objective not supported')
end


%%% EQUALITY CONSTRAINTS %%%
switch optmType
    case {'cmplxAmpl','phase_op1'}
        %keep far-field power constant
        Heq(:,:,1) = 2*eye(numel(x));
    case 'phase_op2'
        Heq = Ha;
    case 'phase_op3'
        Heq=[];
end
%string specifying all other constraints
conStr = extractAfter(Case,'_st:');
if contains(conStr,'Ky=Kx');    Heq = cat(3,Heq, H_Dy-H_Dx );    end
if contains(conStr,'Kxy=Kx');   Heq = cat(3,Heq, H_Dxy-H_Dx );   end
if contains(conStr,'Kxny=Kx');  Heq = cat(3,Heq, H_Dxny-H_Dx );  end
if contains(conStr,'Kxnz=Kxz'); Heq = cat(3,Heq, H_Dxnz-H_Dxz ); end
if contains(conStr,'Kynz=Kyz'); Heq = cat(3,Heq, H_Dynz-H_Dyz ); end
if contains(conStr,'Kynz=Kxz'); Heq = cat(3,Heq, H_Dynz-H_Dxz ); end
if contains(conStr,'Kyz=Kxz');  Heq = cat(3,Heq, H_Dyz-H_Dxz );  end
if contains(conStr,'Kz=Kx*R');  Heq = cat(3,Heq, H_Dz-H_Dx*R );  end
if contains(conStr,'Fx=0');     Heq = cat(3,Heq, H_Bx );         end
if contains(conStr,'Fy=0');     Heq = cat(3,Heq, H_By );         end
if contains(conStr,'Fz=0');     Heq = cat(3,Heq, H_Bz );         end


%%% INEQUALITY CONSTRAINTS %%%
Hineq = [];
if contains(conStr,'Kx<=0');    Hineq = cat(3,Hineq, H_Dx );     end


%---final Hessian of the Lagrangian
if ~isempty(Hineq) 
    LagMultIneq = lambda.ineqnonlin;
    LagMultIneq = reshape(LagMultIneq,1,1,[]);
    LHineq = sum(LagMultIneq.*Hineq,3);
else
    LHineq = zeros(size(Hobj));
end

if ~isempty(Heq) 
    LagMultEq = lambda.eqnonlin;
    LagMultEq = reshape(LagMultEq,1,1,[]);
    LHeq = sum(LagMultEq.*Heq,3);
else
    LHeq = zeros(size(Hobj));
end

Hout = Hobj + LHeq + LHineq;


end


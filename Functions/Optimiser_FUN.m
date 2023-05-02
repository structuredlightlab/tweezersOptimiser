function [OptmSol,grad,hess] = Optimiser_FUN(optmType,Case,Start,OptmSetup,m2v,zxRatio)
% This code is an optimiser for finding light fields with the best
% trapping stiffness for a spherical particle

% INPUT
% optmType      what kind of optimiser do you want to run?
%               'cmplxAmpl': optimises both amplitude and phase
%               'phase_op1': optimises phase indirectly; amplitude is kept constant by rescaling x within the objective, constraint, and Hessian functions
%               'phase_op2': optimises phase indirectly; amplitude is kept constant with constraint functions
%               'phase_op3': optimises phase directly
% Case          string specifying objective and constraints
% Start         string specifying the starting point
% OptmSetup     structure containing all the necessary setup parameters
% m2v           transform matrix from the Bessel to VSWF basis
% zxRatio       requested ratio between z and x stiffness


% OUTPUT
% OptmSol       structure containing all the information about the optimiser solution

warning('off','MATLAB:nearlySingularMatrix')

%---import setup data
unpackOptmSetup;
modes = OptmSetup.Modes;

%---precalculate all the operators
switch optmType
    case {'cmplxAmpl','phase_op1','phase_op2'}
        N_var = N_modes*2;
        Bx   = [Qx 1i*Qx; -1i*Qx Qx];         BxT = transpose(Bx);
        By   = [Qy 1i*Qy; -1i*Qy Qy];         ByT = transpose(By);
        Bz   = [Qz 1i*Qz; -1i*Qz Qz];         BzT = transpose(Bz);
        Dx   = [Kx 1i*Kx; -1i*Kx Kx];         DxT = transpose(Dx);
        Dy   = [Ky 1i*Ky; -1i*Ky Ky];         DyT = transpose(Dy);
        Dz   = [Kz 1i*Kz; -1i*Kz Kz];         DzT = transpose(Dz);
        Dxy  = [Kxy 1i*Kxy; -1i*Kxy Kxy];     DxyT = transpose(Dxy);
        Dxny = [Kxny 1i*Kxny; -1i*Kxny Kxny]; DxnyT = transpose(Dxny);
%         Dxz  = [Kxz 1i*Kxz; -1i*Kxz Kxz];     DxzT = transpose(Dxz);
%         Dxnz = [Kxnz 1i*Kxnz; -1i*Kxnz Kxnz]; DxnzT = transpose(Dxnz);
%         Dyz  = [Kyz 1i*Kyz; -1i*Kyz Kyz];     DyzT = transpose(Dyz);
%         Dynz = [Kynz 1i*Kynz; -1i*Kynz Kynz]; DynzT = transpose(Dynz);
    case 'phase_op3'
        N_var = N_modes;
        %all far-field real-amplitudes are equal, and the powers add up to one
        a = ones(N_var,1); a = a/norm(a); aT = transpose(a);
        Bx   = a.*Qx.*aT;   BxT = transpose(Bx);
        By   = a.*Qy.*aT;   ByT = transpose(By);
        Bz   = a.*Qz.*aT;   BzT = transpose(Bz);
        Dx   = a.*Kx.*aT;   DxT = transpose(Dx);
        Dy   = a.*Ky.*aT;   DyT = transpose(Dy);
        Dz   = a.*Kz.*aT;   DzT = transpose(Dz);
        Dxy  = a.*Kxy.*aT;  DxyT = transpose(Dxy);
        Dxny = a.*Kxny.*aT; DxnyT = transpose(Dxny);
%         Dxz  = a.*Kxz.*aT;  DxzT = transpose(Dxz);
%         Dxnz = a.*Kxnz.*aT; DxnzT = transpose(Dxnz);
%         Dyz  = a.*Kyz.*aT;  DyzT = transpose(Dyz);
%         Dynz = a.*Kynz.*aT; DynzT = transpose(Dynz);
    otherwise
        error('scpecified optimiser type not supported')
end

%---precalculate the Hessian of the amplitude constraint function for the phase-only option 2
if  strcmp(optmType,'phase_op2')
    %all far-field real-amplitudes are equal, and the powers add up to one
    a = ones(N_modes,1); a = a/norm(a);
    Ha = zeros(N_var,N_var,N_modes);
    for i = 1:N_modes; Ha(i,i,i) = 2; Ha(i+N_modes,i+N_modes,i) = 2; end
else
    if strcmp(optmType,'phase_op1') || strcmp(optmType,'cmplxAmpl')
        a = [];
    end
    Ha = [];
end

%---assemble the operators needed for a particular set of objective andconstraint functions
if contains(Case,'Fx');   Ops.Bx = Bx;     Ops.BxT = BxT;     end
if contains(Case,'Fy');   Ops.By = By;     Ops.ByT = ByT;     end
if contains(Case,'Fz');   Ops.Bz = Bz;     Ops.BzT = BzT;     end
if contains(Case,'Kx');   Ops.Dx = Dx;     Ops.DxT = DxT;     end
if contains(Case,'Ky');   Ops.Dy = Dy;     Ops.DyT = DyT;     end
if contains(Case,'Kz');   Ops.Dz = Dz;     Ops.DzT = DzT;     end
if contains(Case,'Kxy');  Ops.Dxy = Dxy;   Ops.DxyT = DxyT;   end
if contains(Case,'Kxny'); Ops.Dxny = Dxny; Ops.DxnyT = DxnyT; end
if contains(Case,'Kxz');  Ops.Dxz = Dxz;   Ops.DxzT = DxzT;   end
if contains(Case,'Kxnz'); Ops.Dxnz = Dxnz; Ops.DxnzT = DxnzT; end
if contains(Case,'Kyz');  Ops.Dyz = Dyz;   Ops.DyzT = DyzT;   end
if contains(Case,'Kynz'); Ops.Dynz = Dynz; Ops.DynzT = DynzT; end


%---define an output function to keep track of optimisation variable values at each iteration
%set options for the solver
constrTol = 1e-5;
optimlTol = 1e-5;
options = optimoptions('fmincon','ConstraintTolerance',constrTol,'UseParallel',false,'Display','final',...
                       'Algorithm','interior-point','SpecifyObjectiveGradient',true,...
                       'SpecifyConstraintGradient',true,...
                       'MaxIterations',1e7,'MaxFunctionEvaluations',1e8,...
                       'StepTolerance',1e-12,'OptimalityTolerance',optimlTol,...
                       'PlotFcn',@PlotFnForOptimiser,...
                       'OutputFcn',@OutputFnForOptimiser,'ScaleProblem','none',...
                       'CheckGradients',false,'FiniteDifferenceType','central',...
                       'HessianFcn',@(x,lambda)hessianfcn(x,lambda,Ops,Case,optmType,zxRatio,Ha));
                   
                   
%create a problem structure
problem.options = options;

%set the objective function
problem.objective = @(x)ObjectiveFun(x,Ops,Case,optmType);

%specify non-linear constraints; see function "nonLinConstraints" for details
problem.nonlcon = @(x)nonLinConstraints(x,Ops,Case,optmType,zxRatio,a);

%specify initial values
switch Start
    case '0oamModes'     %this is a Gaussian-like beam
        switch optmType
            case 'cmplxAmpl'
                problem.x0 = zeros(N_var,1);
                problem.x0(modes(2,:)==0) = 1;
                problem.x0 = problem.x0/norm(problem.x0);
            case 'phase_op2'
                problem.x0 = zeros(N_var,1);
                problem.x0(modes(2,:)==0) = max(a)/sqrt(2); 
                problem.x0 = problem.x0/norm(problem.x0);
            case 'phase_op3'
                problem.x0 = zeros(N_var,1);
%                 problem.x0 = ones(N_var,1);
        end
    case 'random'
        switch optmType
            case 'cmplxAmpl'
                x0_ampl = rand(N_modes,1);
                x0_ampl = x0_ampl/norm(x0_ampl) .* exp(1i*2*pi*rand(N_modes,1));
                problem.x0 = [real(x0_ampl);  imag(x0_ampl)];
                problem.x0 = problem.x0/norm(problem.x0);
            case 'phase_op2'
                x0_ampl = a .* exp(1i*2*pi*rand(N_modes,1));
                problem.x0 = [real(x0_ampl);  imag(x0_ampl)];
                problem.x0 = problem.x0/norm(problem.x0);
            case 'phase_op3'
                problem.x0 = 2*pi*rand(N_var,1);
        end        
    otherwise
        error('specified starting point not supported')
end

%bounds on the optimisation variables
switch optmType
    case {'cmplxAmpl','phase_op1'}
        problem.lb = -ones(N_var,1);
        problem.ub = ones(N_var,1);
    case 'phase_op2'
        problem.lb = -ones(N_var,1)*max(a);
        problem.ub = ones(N_var,1)*max(a);
end

%choose a solver
problem.solver = 'fmincon';

%---Run the optimiser
%solve the problem
[x,feval,exitflag,output,lambda,grad,hess] = fmincon(problem);


%---save the solution
%far-field amplitudes
switch optmType
    case {'cmplxAmpl','phase_op2'}
        sol_E = x(1:end/2) + 1i*x(end/2+1:end);
    case 'phase_op1'
        normF = sqrt(x(1:end/2).^2 + x(end/2+1:end).^2);
        x = x./[normF; normF];
        x = x/sqrt(N_modes);
        sol_E = x(1:end/2) + 1i*x(end/2+1:end);
    case 'phase_op3'
        sol_E = a.*exp(1i*x);
end

%solution in the near-field
sol_bsc = m2v*sol_E;
sol_Beam = ott.Bsc(sol_bsc(1:end/2),sol_bsc(end/2+1:end),'regular','incident');
sol_Beam.wavelength = wavelength0/n_medium;

%save
OptmSol.Case = Case;
OptmSol.StartPoint = Start;
OptmSol.Output = output;
OptmSol.ExitFlag = exitflag;
OptmSol.SolutionInFarField = sol_E;
OptmSol.SolutionInBSCs = sol_Beam;
OptmSol.kxkzRatio = 1/zxRatio;

end


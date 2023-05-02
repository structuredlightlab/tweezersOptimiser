function [H] = hessianGWS(op,opT,x,xConj,xT,optmType)

switch optmType
    case {'cmplxAmpl','phase_op1','phase_op2'}
        H = real(op + opT);
    case 'phase_op3'
        H1 = -diag(xConj).*(op*x);
        H2 = -diag(x).*(opT*xConj);
        H3 = xConj.*op.*xT;
        H = real(H1 + H2 + H3 + H3.');
end


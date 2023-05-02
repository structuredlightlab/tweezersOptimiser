function [Q,K] = calcQK(T,b2v,direction,beam,derivativeType,method)

%INPUT
%T                  T-matrix of the particle
%b2v                transform matrix from Bessel to VSWF basis
%direction          direction along which the GWS operators are evaluated,
%                   specified as 3-element column vector
%beam               a sample Bessel beam
%derivativeType     method for calculating the derivative of the scattering
%                   matrix S; can be:
%                   'finiteDiff3' - central finite difference scheme (uses 2 steps)
%                   'finiteDiff5' - uses 4 steps
%                   'finiteDiff6' - uses 6 steps 
%                   'analytic'    - analytic derivative (valid only for small displacements; fastest method)
%method             the method to be used for evaluating S^(-1); options are: 
%                   'cct'    - complex conjugate transpose (use this method! the others can be numerically unstable) 
%                   'effInv' - effective inverse using singular value decomposition
%                   'solver' - invert directly

%OUTPUT
%Q                  GWS force operator
%K                  stiffness operator



if strcmp(method,'cct'); numSV=0; end

switch derivativeType
    case 'finiteDiff3'
        %obtain translation matrices; these matrices translate the ORIGIN with
        %respect to which the VSWF are defined (and therefore the bead), the light
        %field itself will stay where it was, but its BSCs will be in the new basis
        [~, AB] = translateXyz(beam,direction);      %to +dx 
        [~, CD] = translateXyz(beam,-direction);     %to -dx
        
        %calculate translated T-matrices
        T_pos = T; T_neg = T; 
        T_pos.data = CD*(T_pos.data*AB);
        T_neg.data = AB*(T_neg.data*CD);
        
        %apply the T-matrices to b2v to obtain the scattering matrices
        S_neg  = T_neg.data*b2v;
        S_0    = T.data*b2v;
        S_pos  = T_pos.data*b2v;
        
        dS = (S_pos - S_neg)/(norm(direction)*2);
        d2S = (S_pos - 2*S_0 + S_neg)/(norm(direction)^2);
        S_inv = S_0';
        dS_inv = dS';
        
        Q = -1i*S_inv*dS;
        K = -1i*(dS_inv*dS + S_inv*d2S);
        
    case 'finiteDiff5'
        %obtain translation matrices; these matrices translate the ORIGIN with
        %respect to which the VSWF are defined (and therefore the bead), the light
        %field itself will stay where it was, but its BSCs will be in the new basis
        [~, AB]  = translateXyz(beam,direction);      %to +dx 
        [~, AB2] = translateXyz(beam,direction*2);    %to +2dx
        [~, CD]  = translateXyz(beam,-direction);     %to -dx
        [~, CD2] = translateXyz(beam,-direction*2);   %to -2dx
        
        %calculate translated T_matrices
        T_pos = T; T_neg = T; T_pos2 = T; T_neg2 = T;
        T_pos1.data = CD*(T_pos.data*AB);
        T_neg1.data = AB*(T_neg.data*CD);
        T_pos2.data = CD2*(T_pos2.data*AB2);
        T_neg2.data = AB2*(T_neg2.data*CD2);
                
        %apply the T-matrices to b2v to obtain the scattering matrices
        S_neg2 = T_neg2.data*b2v;
        S_neg  = T_neg1.data*b2v;
        S_0    = T.data*b2v;
        S_pos  = T_pos1.data*b2v;
        S_pos2 = T_pos2.data*b2v;   
        
        %calculate the GWS operators for a given numSV
        Q_pos = GWS(S_pos,S_pos2,S_0,norm(direction),numSV,method);
        Q_neg = GWS(S_neg,S_0,S_neg2,norm(direction),numSV,method);
        %the stiffness operator
        K = (Q_pos-Q_neg)/(norm(direction)*2);
        
    case 'finiteDiff6'
        %obtain translation matrices; these matrices translate the ORIGIN with
        %respect to which the VSWF are defined (and therefore the bead), the light
        %field itself will stay where it was, but its BSCs will be in the new basis
        [~, AB1] = translateXyz(beam,direction/2);      %to +dx/2 
        [~, AB2] = translateXyz(beam,direction);        %to +dx
        [~, AB3] = translateXyz(beam,direction*3/2);    %to +3dx/2
        [~, CD1] = translateXyz(beam,-direction/2);     %to -dx/2
        [~, CD2] = translateXyz(beam,-direction);       %to -dx
        [~, CD3] = translateXyz(beam,-direction*3/2);   %to -3dx/2

        %calculate translated T_matrices
        T_pos1 = T; T_neg1 = T; T_pos2 = T; T_neg2 = T; T_pos3 = T; T_neg3 = T;
        T_pos1.data = CD1*(T_pos1.data*AB1);
        T_neg1.data = AB1*(T_neg1.data*CD1);
        T_pos2.data = CD2*(T_pos2.data*AB2);
        T_neg2.data = AB2*(T_neg2.data*CD2);
        T_pos3.data = CD3*(T_pos3.data*AB3);
        T_neg3.data = AB3*(T_neg3.data*CD3);

        %apply the T-matrices to b2v to obtain the scattering matrices
        S_neg3 = T_neg3.data*b2v;
        S_neg2 = T_neg2.data*b2v;
        S_neg1 = T_neg1.data*b2v;
        S_0    = T.data*b2v;
        S_pos1 = T_pos1.data*b2v;
        S_pos2 = T_pos2.data*b2v;
        S_pos3 = T_pos3.data*b2v;

        %calculate the GWS operators
        Q = GWS(S_0,S_pos2,S_neg2,norm(direction)*2,numSV,method);
        Q_pos = GWS(S_pos2,S_pos3,S_pos1,norm(direction),numSV,method);
        Q_neg = GWS(S_neg2,S_neg1,S_neg3,norm(direction),numSV,method);

        %the stiffness operator
        K = (Q_pos-Q_neg)/(norm(direction)*2);

    case 'analytic'
        [~, dAB_plus] = translateXyz_derivative(beam, 1, direction);
        [~, d2AB_plus] = translateXyz_derivative(beam, 2, direction);

        dT = T.data*dAB_plus - dAB_plus*T.data;
        d2T = 2*(-dAB_plus*T.data*dAB_plus + 0.5*(d2AB_plus*T.data + T.data*d2AB_plus));
        
        S = T.data*b2v;
        dS = dT*b2v;
        d2S = d2T*b2v;
        S_inv = S';
        dS_inv = dS';        

        Q = -1i*S_inv*dS;
        K = -1i*(dS_inv*dS + S_inv*d2S);

    otherwise
        error('case not specified')
end






end


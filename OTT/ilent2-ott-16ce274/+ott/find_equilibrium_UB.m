function [zEq,tr1,tr2] = find_equilibrium_UB(z, fz, findTrapRange, bound)
%FIND_EQUILIBRIUM estimates equilibrium positions from position-force data
%
% zeq = find_equilibrium(z, fz) finds the axial equilibrium given two vectors
% z and fz with the position and force values respectively.
% 
% this function is based on "find_equilibrium"
%
% Author: Une Butaite


% Make sure the vectors are both colum vectors
fz = fz(:);
z = z(:);

signFz = sign(fz);
%find the stable roots of fz (i.e. is the force crossing from +ve to -ve in between neighbouring points?)
isStable = signFz(1:end-1)>signFz(2:end);

%if there are no stable points
if sum(isStable)==0
    zEq = [];
else
    zEq = z(isStable);
    zEqInd = find(isStable);
    %determine the force gradient at the equilibria 
    fzGrad = fz(zEqInd+1) - fz(zEqInd);
    %select the most -ve gradient, i.e. the stiffest one
%     zEq = zEq(fzGrad==min(fzGrad)); 
    %select the equilibrium closest to zero
%     zEq_temp = abs(zEq); zEq = zEq(zEq_temp==min(zEq_temp)); zEq = zEq(1);
%     zEqInd = find(z==zEq);
    %select the equilibrium closest to the middle of the search range
    zEq_temp = abs(zEq-z(end/2)); zEq = zEq(zEq_temp==min(zEq_temp)); zEq = zEq(1);
    zEqInd = find(z==zEq);
    
    %fit a line to find the true equilibrium value
    range = zEqInd:zEqInd+1;
    lineZ = polyfit(z(range),fz(range),1);
    zEq = roots(lineZ);
end

if findTrapRange && ~isempty(zEq)
    %find the points where the force crosses from -ve to +ve
    isUnstable = signFz(1:end-1)<signFz(2:end);

    %find the unstable points to the left of zEq if there are any
    ind = find(isUnstable(1:zEqInd)==1);
    if isempty(ind)
        %if there are no unstable points, consider the trap bound as being
        %infinite
%         tr1 = -inf;
        tr1 = min(z);
    else
        %otherwise identify the unstable point closest to zEq
        ind = ind(end);
        switch bound
            case 'forcePeak'
                %find the location of max force in between zEq and the trap bound
                f = max(fz(ind:zEqInd));
                %find the location of max force in between zEq and the start of range
%                 f = max(fz(1:zEqInd));
                tr1 = z(fz==f);
            case 'forceRoot'
                %fit a line to find the true root
                range = ind:ind+1;
                lineZ = polyfit(z(range),fz(range),1);
                tr1 = roots(lineZ);
%                 tr1 = z(ind);
        end
    end
    
    %find the unstable points to the right of zEq if there are any
    ind = find(isUnstable(zEqInd+1:end)==1);
    if isempty(ind)
        %if there are no unstable points, consider the trap bound as being
        %infinite
%         tr2 = inf;
        tr2 = max(z);
    else
        %otherwise identify the unstable point closest to zEq
        ind = ind(1) + numel(isUnstable(1:zEqInd));
        switch bound
            case 'forcePeak'
                %find the location of min force in between zEq and the trap bound
                f = min(fz(zEqInd+1:ind));
                %find the location of max force in between zEq and the end of range
%                 f = min(fz(zEqInd:end));
                tr2 = z(fz==f);
            case 'forceRoot'
                %fit a line to find the true root
                range = ind:ind+1;
                lineZ = polyfit(z(range),fz(range),1);
                tr2 = roots(lineZ);
%                 tr2 = z(ind);
        end
    end
    
else
    tr1=[]; tr2=[];
end

end

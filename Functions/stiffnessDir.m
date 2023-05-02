function [kappa] = stiffnessDir(directions,beam,T,eq)
%stiffnessDir calculates the optical trapping stiffness along specified
%directions

%directions         each column is the direction (which takes into account
%                   the step size) along which stiffness is estimated, m
%beam               ott.Bsc structure describing the trapping beam
%T                  T-matrix of the trapped particle
%eq                 column vector with the equilibrium location, m


kappa = zeros(1,size(directions,2));
for i = 1:size(directions,2)
    dir = directions(:,i);
    dirMag = norm(dir);
    
    %force components along x,y,z; 1st column is on the +ve side of eq, 2nd - on the -ve side
    forceXYZ = ott.forcetorque(beam,T,'position',[dir, -dir]+eq);
    %force magnitude along the specified direction
    forceDir = dot(forceXYZ,[1 1].*dir/dirMag);
    
    kappa(i) = (forceDir(1)-forceDir(2))/(2*dirMag);
end

end


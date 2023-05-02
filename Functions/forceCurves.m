function [force,pos] = forceCurves(dir,range,Np,beams,eq,T)


%no. of beams
Nb = numel(beams);
%no. of directions
Nd = size(dir,2);

%points at which force is evaluated
pos = linspace(-1,1,Np)*range;
force = cell(Nd,Nb);
for i = 1:Nb
   for j = 1:Nd
       d = dir(:,j);
       dirMag = norm(d);
       f = ott.forcetorque(beams(i), T, 'position', d*pos+eq(:,i));
       force{j,i} = dot(f,ones(1,Np).*d/dirMag);
   end
end


end


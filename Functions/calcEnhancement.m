function [enh_o,enh_eq,k,k_ref_o,k_ref_eq,zEq,forceGz] = calcEnhancement(beam,beamRef,directions,dx,radius,T,atOrgn,atEquil)

%directions         axes along which to estimate enhancment; each column
%                   contains a three-element vector defining  the direction in Cartesian
%                   coordinates, and does not need to be normalised
%radius             radius of the trapped particle, m
%T                  T-matrix
%atOrgn             if = true, calculates enhancement at the origin
%atEquil            if = true, calculates enhancement at the equilibrium 
%                   (for reference beam equilibrium might be different from
%                   the origin)

%enhancement        stiffness enhancement factor along every direction. If
%                   no equilibrium exists, the function returns Nan


%default values
enh_o = zeros(1,size(directions,2))*nan;
enh_eq = zeros(1,size(directions,2))*nan;
k_ref_o = zeros(1,size(directions,2))*nan;
k_ref_eq = zeros(1,size(directions,2))*nan;

%normalise all directions and scale with the step size
dirs = directions./vecnorm(directions)*dx;

%equilibrium of the optimised beam
eq = [0;0;0];
%stiffness in every direction for the original beam
k = stiffnessDir(dirs,beam,T,eq);

%-----stiffness of the reference beam
%at the origin
if atOrgn
    eq_ref = [0;0;0];
    k_ref_o = stiffnessDir(dirs,beamRef,T,eq_ref);
    %calculate stiffness enhancment
    enh_o = k./k_ref_o;
end
%at the equilibrium
if atEquil
    %find the z-equilibrium of the reference beam
    Np = 200;     rng = linspace(-radius*2,radius*2,Np);
    forceGz = ott.forcetorque(beamRef,T,'position',[0;0;1]*rng);
    zEq = ott.find_equilibrium_UB(rng,forceGz(3,:),0,'forcePeak');
    eq_ref = [0;0;zEq];
    
    if ~isempty(zEq)
        k_ref_eq = stiffnessDir(dirs,beamRef,T,eq_ref);
        %calculate stiffness enhancment
        enh_eq = k./k_ref_eq;
    else
        zEq = nan;
    end
    
end







end


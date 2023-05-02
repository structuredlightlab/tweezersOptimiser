function [I_xy,I_xz] = plotIntCrossSections(beam_optim,T_sc,T_int,r_plot,res,radius)


%scattered beam
beam_sc = T_sc*beam_optim;

%total outgoing beam
beam_out = beam_sc.totalField(beam_optim);
beam_out.basis = 'regular';

%beam inside the particle
beam_int =  T_int*beam_optim;



%-----XY cross-section-----
%coordinates
range = linspace(-r_plot, r_plot, res);
[x, y, z] = meshgrid(range, range, 0);
xy = [x(:) y(:) z(:)].';
%indices where the particle is
idx_particle = (xy(1,:).^2 + xy(2,:).^2 + xy(3,:).^2) < radius^2;

%the complex field
E_xy = zeros(size(xy));
E_xy(:,~idx_particle) = beam_out.emFieldXyz(xy(:,~idx_particle));
E_xy(:,idx_particle)  = beam_int.emFieldXyz(xy(:,idx_particle));

%intensity
I_xy = sum(abs(E_xy).^2,1);
I_xy = reshape(I_xy,[res,res,1]);



%-----XZ cross-section-----
%coordinates
[x, y, z] = meshgrid(range, 0, range);
xz = [x(:) y(:) z(:)].';
%indices where the particle is
idx_particle = (xz(1,:).^2 + xz(2,:).^2 + xz(3,:).^2) < radius^2;

%the complex field
E_xz = zeros(size(xz));
E_xz(:,~idx_particle) = beam_out.emFieldXyz(xz(:,~idx_particle));
E_xz(:,idx_particle)  = beam_int.emFieldXyz(xz(:,idx_particle));

%intensity
I_xz = sum(abs(E_xz).^2,1);
I_xz = reshape(I_xz,[res,res,1]).';


%-----figure-----
%plot the intensitiy of the solution field in the two planes
figure('Position',[200 200 420 220]);
ha = tight_subplot(1,2,.02,[.01 .1],.01);
axes(ha(1))
    imagesc(range, range, I_xy)
    markBead(radius)
    title('Intensity; xy cross-section'); 
    axis equal tight off
    set(gca,'YDir','normal'); 
axes(ha(2))
    imagesc(range, range, I_xz)
    markBead(radius)
    title('Intensity; xz cross-section'); 
    axis equal tight off
    set(gca,'YDir','normal'); 
colormap('gray')



function [] = markBead(radius)
    pos = radius*[-1 -1 2 2];
    rectangle('position',pos,'curvature',[1 1],'edgecolor','w')
    rectangle('position',pos,'curvature',[1 1],'edgecolor','r','linestyle','--')
end

end


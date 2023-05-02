function [refGauss, E] = referenceGaussian(N_modes,modes,b2v,wavelength_medium)
%this function generates a reference Gaussian beam (made up from Bessel beams) 
%for calculating enhancement

%far-field amplitudes
E = zeros(N_modes,1);
E(modes(2,:)==0) = 1;
% theta = unique(modes(1,:));  E(modes(2,:)==0) = exp(-sin(theta).^2);
E = E/norm(E);

bsc = b2v*E;
refGauss = ott.Bsc(bsc(1:end/2),bsc(end/2+1:end),'regular','incident');
refGauss.wavelength = wavelength_medium;



end


function [bsc] = BesselFiniteWidth(alpha,dAlpha,Nb,L,NB,nmax,polrz)

%INPUT
%alpha                 cone angle of the Bessel beam, rad
%dAlpha                angular width of the Bessel beam, rad
%Nb                    number of infinitesimally thin Bessels making up the main Bessel
%L                     OAM of the Bessel
%NB                    total number of modes in the far field
%nmax                  cut-off value for the VSWF expansion
%polrz                 polarisation vector, 2-element vector

%OUTPUT
%bsc                    beam shape coefficients of the finite width Bessel beam; 1D array

%make sure that Nb is odd (in order to use the composite Simpson's rule)
if rem(Nb,2)==0
    Nb = Nb+1; 
    warning('N_b should be odd, added 1 to N_b to correct') 
end

bsc = zeros(sum(3:2:nmax*2+1)*2,Nb);
%cone angles of the infinitesimally thin Bessels
alphaSub = linspace(alpha-dAlpha,alpha,Nb);
%spacing between the infinitesimally thin Bessels
dalpha = alphaSub(2)-alphaSub(1);
for j=1:Nb
        %generate an infinitesimally thin Bessel beam
        besselbeam = ott.BscBessel(nmax, alphaSub(j),'polarisation', polrz,'lmode', L);
        %normalise
        besselbeam = besselbeam*sqrt(sin(alphaSub(j)));
        bsc(:,j) = [besselbeam.a; besselbeam.b];
end
%composite Simpson's rule for the integral
SimpsCoeff = 1/3*[1 4 repmat([2,4],1,(Nb(1)-3)/2) 1];
bsc = bsc.*SimpsCoeff;     
%final normalisation
% bsc = sum(bsc,2)*sqrt(dalpha/(NB*Nb));
bsc = sum(bsc,2)*sqrt(dalpha/(Nb));         %don't need to divide by NB if the far field is already normalised

end


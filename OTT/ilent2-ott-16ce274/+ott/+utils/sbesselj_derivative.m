function [jn_derivative] = sbesselj_derivative(n,kr,derivativeOrder)
% SBESSELJ_DERIVATIVE calculates the 1st or 2nd order derivatives of a
% spherical bessel function jn(kr) in the limit kr->0
% jn(kr) = sqrt(pi/2kr) Jn+0.5(kr)
% 
% This function is NOT part of the optical tweezers toolbox.
% Author: Une B

import ott.*
ott.warning('internal');

kr=kr(:);
n=n(:);

[n,kr]=meshgrid(n,kr);
jn_derivative = zeros(size(n));

if derivativeOrder==1
    jn_derivative(n==1) = 1/3;
elseif derivativeOrder==2
    jn_derivative(n==0) = -1/3;
    jn_derivative(n==2) = 2/15;
else
    error('Only derivatives up to 2nd order are supported')
end

ott.warning('external');

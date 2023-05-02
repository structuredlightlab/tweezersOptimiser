function [A,B,C,Cof0] = translate_z_derivative(nmax,z, derivativeOrder, varargin)
% Calculates derivatives of translation matrices for translation of VSWFs along z axis.
%
% Parameters
%   - nmax (int) -- Determines the number of multipole terms to include
%     in the translation matrices (multipole order).  Can be a single
%     integer or two integers for the ``[row, column]`` nmax.
%     If the row, column indices don't match, A and B will not be square.
%   - z (numeric) -- Translation distance.
%
% This function is NOT part of the optical tweezers toolbox, but it was adapted from OTT's "translate_z".
% Author: Une B


import ott.utils.*
ott.warning('internal');

p = inputParser;
p.addParameter('type', 'sbesselj');
p.addParameter('method', 'gumerov');
p.parse(varargin{:});

if numel(z)>1
    A=cell(numel(z),1);
    B=A;
    for ii=1:numel(z)
        [A{ii},B{ii}]=translate_z_derivative(nmax,z(ii), varargin{:});
    end
    C=0;
    ott.warning('external');
    return
end

% Calculate Nmax for each dimension
if numel(nmax) == 1
  nmax1 = nmax;
  nmax2 = nmax;
else
  nmax1 = nmax(1);
  nmax2 = nmax(2);
  nmax = max(nmax(1:2));
end

if z==0
    A=speye(nmax1^2+nmax1*2, nmax2^2+nmax2*2);
    B=sparse(nmax1^2+nmax1*2, nmax2^2+nmax2*2);
    C=A;
    ott.warning('external');
    return
end

% Calculate the scalar coefficients
switch derivativeOrder
    case 1
        Cof0 = translate_z_gumerov_0(nmax1, nmax2, nmax, 1e-50, p);
        dC = 2*pi*   translate_z_gumerov_derivative(nmax1, nmax2, nmax, abs(z), p, 1);
        [A,B] = calculate_AB_1stDerivative(Cof0, dC, nmax1, nmax2, nmax, z, p);
        C = dC;
    case 2
        dC = 2*pi*   translate_z_gumerov_derivative(nmax1, nmax2, nmax, abs(z), p, 1);
        d2C = (2*pi)^2*translate_z_gumerov_derivative(nmax1, nmax2, nmax, abs(z), p, 2);
        [A,B] = calculate_AB_2ndDerivative(dC, d2C, nmax1, nmax2, nmax, z, p);
        C = d2C;
    otherwise
        error('only derivatives up to 2nd order are supported')
end

if nargout>2
    C=C(1:nmax+1,1:nmax+1,1:min(nmax1, nmax2)+1);
end

ott.warning('external');

end

function [A,B] = calculate_AB_1stDerivative(Cof0, dC, nmax1, nmax2, nmax, z, p)

import ott.utils.*

% OK, that's the scalar coefficients
% Time to find the vector coefficients - Videen (43) & (44)

nn = 1:nmax1;
kk = (1:nmax2).';

matrixm=sqrt(kk.*(kk+1)) ./ sqrt(nn.*(nn+1));

central_iterator1=[1:nmax1].*[2:nmax1+1];
central_iterator2=[1:nmax2].*[2:nmax2+1];

[ciy,cix]=meshgrid(central_iterator1,central_iterator2);

mmm=0;

Cp = Cof0(3:(nmax2+2),2:(nmax1+1),mmm+1);
Cm = Cof0(1:nmax2,2:(nmax1+1),mmm+1);

dC0 = dC(2:(nmax2+1),2:(nmax1+1),mmm+1);

t = matrixm.*(dC0 - 2*pi./(kk+1) .* ...
    sqrt((kk-mmm+1).*(kk+mmm+1)./((2*kk+1).*(2*kk+3))) .* Cp - ...
    2*pi./kk.*sqrt((kk-mmm).*(kk+mmm)./((2*kk+1).*(2*kk-1))).*Cm);

toIndexy=(ciy(:));
toIndexx=(cix(:));
A=t(:);
B=zeros(size(A));

% Total size of A and B: sum((1:nmax).^2)*2 + nmax^2

for mmm=1:min(nmax1, nmax2)

    sz1 = mmm:nmax2;
    sz2 = mmm:nmax1;

    C0 = Cof0((1+mmm):(nmax2+1),(1+mmm):(nmax1+1),mmm+1);
    Cp = Cof0((2+mmm):(nmax2+2),(1+mmm):(nmax1+1),mmm+1);
    Cm = Cof0((mmm):nmax2,(1+mmm):(nmax1+1),mmm+1);
    
    dC0 = dC((1+mmm):(nmax2+1),(1+mmm):(nmax1+1),mmm+1);

    tt = matrixm(sz1, sz2).*(dC0 - 2*pi./(kk(sz1)+1) .* ...
        sqrt((kk(sz1)-mmm+1).*(kk(sz1)+mmm+1) ...
        ./((2*kk(sz1)+1).*(2*kk(sz1)+3))) .* Cp - ...
        2*pi./kk(sz1).*sqrt((kk(sz1)-mmm) ...
        .*(kk(sz1)+mmm)./((2*kk(sz1)+1).*(2*kk(sz1)-1))).*Cm);

    ciys=ciy(mmm:end,mmm:end);
    cixs=cix(mmm:end,mmm:end);

    toIndexy=[toIndexy;(ciys(:)+mmm);(ciys(:)-mmm)];
    toIndexx=[toIndexx;(cixs(:)+mmm);(cixs(:)-mmm)];
    A=[A;tt(:);tt(:)];

    tt = mmm./(kk(sz1).*(kk(sz1)+1)).*matrixm(sz1, sz2) .* C0;
    B=[B;tt(:);-tt(:)];

end

% Keep B real until the end, makes things run faster
B = 1i*2*pi*B;

% This is faster than A = A + sparse(...) and A(sub2ind(...)) = [...]
if z < 0
  [n1, ~] = ott.utils.combined_index(toIndexy);
  [n2, ~] = ott.utils.combined_index(toIndexx);
  B=sparse(toIndexy,toIndexx,B.*(-1).^(n1-n2+1),nmax1*(nmax1+2),nmax2*(nmax2+2));
  A=sparse(toIndexy,toIndexx,A.*(-1).^(n1-n2),nmax1*(nmax1+2),nmax2*(nmax2+2));
else
  B=sparse(toIndexy,toIndexx,B,nmax1*(nmax1+2),nmax2*(nmax2+2));
  A=sparse(toIndexy,toIndexx,A,nmax1*(nmax1+2),nmax2*(nmax2+2));
end

end % calculate_AB

function [A,B] = calculate_AB_2ndDerivative(dC, d2C, nmax1, nmax2, nmax, z, p)

import ott.utils.*

% OK, that's the scalar coefficients
% Time to find the vector coefficients - Videen (43) & (44)

nn = 1:nmax1;
kk = (1:nmax2).';

matrixm=sqrt(kk.*(kk+1)) ./ sqrt(nn.*(nn+1));

central_iterator1=[1:nmax1].*[2:nmax1+1];
central_iterator2=[1:nmax2].*[2:nmax2+1];

[ciy,cix]=meshgrid(central_iterator1,central_iterator2);

mmm=0;

dCp = dC(3:(nmax2+2),2:(nmax1+1),mmm+1);
dCm = dC(1:nmax2,2:(nmax1+1),mmm+1);

d2C0 = d2C(2:(nmax2+1),2:(nmax1+1),mmm+1);

t = matrixm.*(d2C0 - 2*pi./(kk+1) .* ...
    sqrt((kk-mmm+1).*(kk+mmm+1)./((2*kk+1).*(2*kk+3))) .* dCp*2 - ...
    2*pi./kk.*sqrt((kk-mmm).*(kk+mmm)./((2*kk+1).*(2*kk-1))).*dCm*2);

toIndexy=(ciy(:));
toIndexx=(cix(:));
A=t(:);
B=zeros(size(A));

% Total size of A and B: sum((1:nmax).^2)*2 + nmax^2

for mmm=1:min(nmax1, nmax2)

    sz1 = mmm:nmax2;
    sz2 = mmm:nmax1;
    
    dC0 = dC((1+mmm):(nmax2+1),(1+mmm):(nmax1+1),mmm+1);
    dCp = dC((2+mmm):(nmax2+2),(1+mmm):(nmax1+1),mmm+1);
    dCm = dC((mmm):nmax2,(1+mmm):(nmax1+1),mmm+1);
    
    d2C0 = d2C((1+mmm):(nmax2+1),(1+mmm):(nmax1+1),mmm+1);
    
    tt = matrixm(sz1, sz2).*(d2C0 - 2*pi./(kk(sz1)+1) .* ...
        sqrt((kk(sz1)-mmm+1).*(kk(sz1)+mmm+1) ...
        ./((2*kk(sz1)+1).*(2*kk(sz1)+3))) .* dCp*2 - ...
        2*pi./kk(sz1).*sqrt((kk(sz1)-mmm) ...
        .*(kk(sz1)+mmm)./((2*kk(sz1)+1).*(2*kk(sz1)-1))).*dCm*2);

    ciys=ciy(mmm:end,mmm:end);
    cixs=cix(mmm:end,mmm:end);

    toIndexy=[toIndexy;(ciys(:)+mmm);(ciys(:)-mmm)];
    toIndexx=[toIndexx;(cixs(:)+mmm);(cixs(:)-mmm)];
    A=[A;tt(:);tt(:)];

    tt = mmm./(kk(sz1).*(kk(sz1)+1)).*matrixm(sz1, sz2) .* dC0;
    B=[B;tt(:);-tt(:)];

end

% Keep B real until the end, makes things run faster
B = 1i*2*pi*B*2;

% This is faster than A = A + sparse(...) and A(sub2ind(...)) = [...]
if z < 0
  [n1, ~] = ott.utils.combined_index(toIndexy);
  [n2, ~] = ott.utils.combined_index(toIndexx);
  B=sparse(toIndexy,toIndexx,B.*(-1).^(n1-n2+1),nmax1*(nmax1+2),nmax2*(nmax2+2));
  A=sparse(toIndexy,toIndexx,A.*(-1).^(n1-n2),nmax1*(nmax1+2),nmax2*(nmax2+2));
else
  B=sparse(toIndexy,toIndexx,B,nmax1*(nmax1+2),nmax2*(nmax2+2));
  A=sparse(toIndexy,toIndexx,A,nmax1*(nmax1+2),nmax2*(nmax2+2));
end

end % calculate_AB

function C = translate_z_gumerov_derivative(nmax1, nmax2, nmax, r, p, derivativeOrder)
% Optimised implementation from Alex

import ott.utils.*

% Having pre-computed a_nm it's fast?
% nmax=3;
% r=0.1;

%some "setup" values:
mmax=min(nmax1,nmax2);
m=0;
fval=2*nmax+1;
nd=[m:fval];

kr=2*pi*r;

%compute seed functions:
switch p.Results.type
  case 'sbesselj'       % regular to regular
    C_nd00=[sqrt(2*nd+1).*ott.utils.sbesselj_derivative(nd,kr,derivativeOrder)];
  otherwise
      error('derivatives only supported for regular to regular translations')  
end

C_ndn0=zeros(length(nd)+1,length(nd)+1);
C_ndn0(1+[1:length(C_nd00)],2)=C_nd00;
C_ndn0(2,1+[1:length(C_nd00)])=((-1).^(nd).*C_nd00).';

%gumerov's zonal coefficients are m=0. Compute columns, limited by diagonal:
%compute lower diagonal first:
for jj=1:nmax
    ii=[jj:fval-jj].';
    C_ndn0(ii+2,ii(1)+2)=(anm_l(ii(1)-2,0).*C_ndn0(ii+2,ii(1))-anm_l(ii,0).*C_ndn0(ii+3,ii(1)+1)+anm_l(ii-1,0).*C_ndn0(ii+1,ii(1)+1))./anm_l(ii(1)-1,0);
    C_ndn0(ii(1)+2,ii+2)=((-1).^(jj+ii).*C_ndn0(ii+2,ii(1)+2)).';
end

%create "C":
C=zeros(nmax2+2,nmax1+1,mmax+1);
C(:,:,1)=C_ndn0(2:(nmax2+3),2:(nmax1+2));

%Having computed anm for m=0; cases we now can compute anm for all
%remaining cases:
ANM=anm_l([0:2*nmax+1].',[1:nmax]);
IANM=1./ANM;
for m=1:mmax

    %having computed the zonal coefficients we now compute the "diagonal ones"
    %(tesseral)
    %i.e. ones which generate m on the first column we then reproduce the same
    %commputation for the n nd recursion:

    nd=[m:fval-m].';
    C_nd1m=(bnm_l(nd,-m).*C_ndn0(nd+1,m+1)-bnm_l(nd+1,m-1).*C_ndn0(nd+3,m+1))./bnm_l(m,(-m));

    %having computed the first seed column we now recur the elements:
    C_ndn1=zeros(size(C_ndn0)); %make zero as we re-use
    C_ndn1([1:length(C_nd1m)]+m+1,m+2)=C_nd1m;
    C_ndn1(m+2,[1:length(C_nd1m)]+m+1)=((-1).^(nd+m).*C_nd1m).';

    for jj=m+1:nmax
        ii=[jj:fval-jj].';
%         C_ndn1(ii+2,ii(1)+2)=(anm(ii(1)-2,m).*C_ndn1(ii+2,ii(1))-anm(ii,m).*C_ndn1(ii+3,ii(1)+1)+anm(ii-1,m).*C_ndn1(ii+1,ii(1)+1))./anm(ii(1)-1,m);
        C_ndn1(ii+2,ii(1)+2)=(ANM(ii(1)-1,m).*C_ndn1(ii+2,ii(1))-ANM(ii+1,m).*C_ndn1(ii+3,ii(1)+1)+ANM(ii,m).*C_ndn1(ii+1,ii(1)+1)).*IANM(ii(1),m);
        C_ndn1(ii(1)+2,ii+2)=((-1).^(jj+ii).*C_ndn1(ii+2,ii(1)+2)).';
    end
    C_ndn0=C_ndn1;

    C(:,:,m+1)=C_ndn0(2:(nmax2+3),2:(nmax1+2));

end

end % translate_z_gumerov

function C = translate_z_gumerov_0(nmax1, nmax2, nmax, r, p)
% this is a copy of translate_z_gumerov

% Optimised implementation from Alex

import ott.utils.*

% Having pre-computed a_nm it's fast?
% nmax=3;
% r=0.1;

%some "setup" values:
mmax=min(nmax1,nmax2);
m=0;
fval=2*nmax+1;
nd=[m:fval];

kr=2*pi*r;

%compute seed functions:
switch p.Results.type
  case 'sbesselj'       % regular to regular
    C_nd00=[sqrt(2*nd+1).*ott.utils.sbesselj(nd,kr)];
  otherwise
      error('translate_z_gumerov_0 only supports regular to regular translations')
end

C_ndn0=zeros(length(nd)+1,length(nd)+1);
C_ndn0(1+[1:length(C_nd00)],2)=C_nd00;
C_ndn0(2,1+[1:length(C_nd00)])=((-1).^(nd).*C_nd00).';

%gumerov's zonal coefficients are m=0. Compute columns, limited by diagonal:
%compute lower diagonal first:
for jj=1:nmax
    ii=[jj:fval-jj].';
    C_ndn0(ii+2,ii(1)+2)=(anm_l(ii(1)-2,0).*C_ndn0(ii+2,ii(1))-anm_l(ii,0).*C_ndn0(ii+3,ii(1)+1)+anm_l(ii-1,0).*C_ndn0(ii+1,ii(1)+1))./anm_l(ii(1)-1,0);
    C_ndn0(ii(1)+2,ii+2)=((-1).^(jj+ii).*C_ndn0(ii+2,ii(1)+2)).';
end

%create "C":
C=zeros(nmax2+2,nmax1+1,mmax+1);
C(:,:,1)=C_ndn0(2:(nmax2+3),2:(nmax1+2));

%Having computed anm for m=0; cases we now can compute anm for all
%remaining cases:
ANM=anm_l([0:2*nmax+1].',[1:nmax]);
IANM=1./ANM;
for m=1:mmax

    %having computed the zonal coefficients we now compute the "diagonal ones"
    %(tesseral)
    %i.e. ones which generate m on the first column we then reproduce the same
    %commputation for the n nd recursion:

    nd=[m:fval-m].';
    C_nd1m=(bnm_l(nd,-m).*C_ndn0(nd+1,m+1)-bnm_l(nd+1,m-1).*C_ndn0(nd+3,m+1))./bnm_l(m,(-m));

    %having computed the first seed column we now recur the elements:
    C_ndn1=zeros(size(C_ndn0)); %make zero as we re-use
    C_ndn1([1:length(C_nd1m)]+m+1,m+2)=C_nd1m;
    C_ndn1(m+2,[1:length(C_nd1m)]+m+1)=((-1).^(nd+m).*C_nd1m).';

    for jj=m+1:nmax
        ii=[jj:fval-jj].';
%         C_ndn1(ii+2,ii(1)+2)=(anm(ii(1)-2,m).*C_ndn1(ii+2,ii(1))-anm(ii,m).*C_ndn1(ii+3,ii(1)+1)+anm(ii-1,m).*C_ndn1(ii+1,ii(1)+1))./anm(ii(1)-1,m);
        C_ndn1(ii+2,ii(1)+2)=(ANM(ii(1)-1,m).*C_ndn1(ii+2,ii(1))-ANM(ii+1,m).*C_ndn1(ii+3,ii(1)+1)+ANM(ii,m).*C_ndn1(ii+1,ii(1)+1)).*IANM(ii(1),m);
        C_ndn1(ii(1)+2,ii+2)=((-1).^(jj+ii).*C_ndn1(ii+2,ii(1)+2)).';
    end
    C_ndn0=C_ndn1;

    C(:,:,m+1)=C_ndn0(2:(nmax2+3),2:(nmax1+2));

end

end % translate_z_gumerov

function a_nm = anm_l(n,m)
% For translate_z_gumerov
fn=1./(2*n+1)./(2*n+3);
a_nm=sqrt((n+abs(m)+1).*(n-abs(m)+1).*fn);
a_nm(n<0)=0;
a_nm(abs(m)>n)=0;
end

function b_nm = bnm_l(n,m)
% For translate_z_gumerov
b_nm=(2*(m<0)-1).*sqrt((n-m-1).*(n-m)./(2*n-1)./(2*n+1));
b_nm(abs(m)>n)=0;
end


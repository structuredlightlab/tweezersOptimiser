
InputBasis = OptmSetup.InputBasis;
modes = OptmSetup.Modes;
N_modes = OptmSetup.Nmodes;

if isfield(OptmSetup,'m2v');   m2v = OptmSetup.modes2vswf; end

if isfield(OptmSetup,'Qx');    Qx = OptmSetup.Qx;  end
if isfield(OptmSetup,'Qy');    Qy = OptmSetup.Qy;  end
if isfield(OptmSetup,'Qxy');   Qxy = OptmSetup.Qxy;  end
if isfield(OptmSetup,'Qxny');  Qxny = OptmSetup.Qxny;  end
if isfield(OptmSetup,'Qz');    Qz = OptmSetup.Qz;  end
if isfield(OptmSetup,'Kx');    Kx = OptmSetup.Kx;  end
if isfield(OptmSetup,'Ky');    Ky = OptmSetup.Ky;  end
if isfield(OptmSetup,'Kxy');   Kxy = OptmSetup.Kxy;  end
if isfield(OptmSetup,'Kxny');  Kxny = OptmSetup.Kxny;  end
if isfield(OptmSetup,'Kxz');   Kxz = OptmSetup.Kxz;  end
if isfield(OptmSetup,'Kxnz');  Kxnz = OptmSetup.Kxnz;  end
if isfield(OptmSetup,'Kyz');   Kyz = OptmSetup.Kyz;  end
if isfield(OptmSetup,'Kynz');  Kynz = OptmSetup.Kynz;  end
if isfield(OptmSetup,'Kz');    Kz = OptmSetup.Kz;  end

NA = OptmSetup.NA;
n_medium = OptmSetup.nMedium;
n_particle = OptmSetup.nParticle;
wavelength0 = OptmSetup.wavelength0;
wavelength_medium = wavelength0/n_medium;
polrz = OptmSetup.polarisation;
radius = OptmSetup.particleRadius;
nmax = OptmSetup.nMax;

T = OptmSetup.Tmatrix;
# tweezersOptimiser
These are the codes that were used to simulate the results published in the "Photon-efficient optical tweezers via wavefront shaping" paper.

The "madeToMeasure" code contains everything you need to tailor your own optical trap for a chosen spherical particle. 

The "Functions" folder contains functions called by "madeToMeasure".

The "OTT" folder contains the Optical Tweezers Toolbox codes in the exact version that was used for this project, which includes a few minor modifications and a few custom functions; see below for more details. 
OTT is copyright 2007-2018 by The University of Queensland under the Creative Commons Attribution-NonCommercial 4.0 International License http://creativecommons.org/licenses/by-nc/4.0/. Documentation for OTT is available at https://ott.readthedocs.io/en/latest/. <br/><br/><br/><br/>



These are the functions in the "OTT" folder that are not part of the OTT suite and were written by us:  
find_equilibrium_UB  
sbesselj_derivative  
translate_z_derivative  
translateZ_type_helper_derivative  
translateXyz_derivative  
translateRtp_derivative

The following changes were made to OTT functions:  
in BscPmGauss line 128 changed the default of "angular_scaling" to be 'sintheta' instead of 'tantheta';  
in BscBessel lines 153 & 157 removed the "-" sign in front of i to change the default direction of beam propagation;  
in BscBessel line 170 made sure that "beam" output is never empty.




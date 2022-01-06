% CALC_SHEAR Calculate shear for polar data.
%   [U0,UR,UTHETA,SHMAG,URES] = CALC_SHEAR(R,AZ,V,CONF) computes
%   the "base value" U0, radial shear UR, azimuthal shear UTHETA,
%   shear magnitude SHMAG, and residual field URES for the field
%   data (e.g., velocities) having ranges R, azimuths AZ (degrees), 
%   values V (ranges down rows, azimuths across columns) and 
%   confidence weights CONF.  The shears are computed via a 
%   CONF-weighted least-squares linear fit (see WTDREGRESS).  The 
%   outputs U0 and URES have the units of V, while UR, UTHETA and
%   SHMAG have units (units of V)/(units of R).
%
%   Example:
%
%    r = 5000:100:5500; 
%    az = 100:2:112;
%    [dAz,R] = meshgrid((az-mean(az)).*pi/180,r);
%    noise = randn(size(R));
%    conf = normpdf(noise,0,0.2);
%    v = 3 - 0.004*R.*dAz + 0.008*(R-mean(r)) + noise
%    [u0,ur,utheta,shmag,ures] = calc_shear(r,az,v,conf)
%
%  See also: WTDREGRESS, REGRESS, WTDPOLYFIT.

%  Written by John K. Williams (jkwillia@ucar.edu, 303-497-2822)


function [u0,ur,utheta,shmag,ures] = calc_shear(r,az,v,conf)

  az = az*pi/180;
  [dAz,R] = meshgrid(az-mean(az),r);
  mat = [ones(prod(size(v)),1), R(:)-mean(r), R(:).*dAz(:)];
  b = wtdregress(v(:),mat,conf(:));
  
  u0 = b(1);
  ur = b(2);
  utheta = b(3);
  shmag = sqrt(sum(b(2:3).^2));
  ures = reshape(v(:)-mat*b,size(v));
  
% END (calc_shear)
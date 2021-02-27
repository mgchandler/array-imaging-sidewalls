function scat_amp = fn_point_scatterer_amp(in_mode, out_mode, long_speed, shear_speed)
% Computes the scattering amplitude for a point scatterer.
%
% INPUTS:
% - in_mode : logical
%       Logical test of whether the incident mode on the scatterer is
%       transverse.
% - out_mode : logical
%       Logical test of whether the output mode from the scatterer is
%       transverse.
% - long_speed : double
%       Longitudinal speed in the block material.
% - shear_speed : double
%       Shear speed in the block material.
%
% OUTPUTS:
% - scat_amp : complex
%       Scattering amplitude from a point scatterer. Need to adapt this
%       function to allow matrices to be output to speed up function. Note
%       that this function should be the easiest to make this change!

if and(~in_mode, ~out_mode) % If view is L-L
    scat_amp = 1.0;
elseif and(~in_mode, out_mode) % If view is L-T
    scat_amp = long_speed / shear_speed;
elseif and(in_mode, ~out_mode) % If view is T-L
    scat_amp = - shear_speed / long_speed;
elseif and(in_mode, out_mode) % If view is T-T
    scat_amp = 1.0;
end

end
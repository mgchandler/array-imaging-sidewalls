function val = fn_scattering_bilinear_interp(scat_matrix, inc_angle, out_angle)
% Finds the scattering amplitude from inc_angle to out_angle from the
% scattering matrix provided using bilinear interpolation.
%
% INPUTS:
% - scat_matrix : complex square array
%       Pre-computed scattering matrix providing scattering amplitudes when
%       a ray is scattered from inc_angle to out_angle.
% - inc_angle : double
%       Incident angle in range [-pi, pi).
% - out_angle : double
%       Outgoing angle in range [-pi, pi).
%
% OUTPUTS:
% - val : complex
%       The interpolated value obtained from the matrix.

[num_points, ~] = size(scat_matrix);
dtheta = 2*pi / num_points;

inc_angle_idx = mod(floor((inc_angle + pi) / dtheta), num_points);
out_angle_idx = mod(floor((out_angle + pi) / dtheta), num_points);

inc_angle_frac = mod(inc_angle + pi, dtheta) / dtheta;
out_angle_frac = mod(out_angle + pi, dtheta) / dtheta;

inc_angle_idx_plus1 = mod(inc_angle_idx + 1, num_points);
out_angle_idx_plus1 = mod(out_angle_idx + 1, num_points);

inc_angle_idx = inc_angle_idx + 1;
out_angle_idx = out_angle_idx + 1;
inc_angle_idx_plus1 = inc_angle_idx_plus1 + 1;
out_angle_idx_plus1 = out_angle_idx_plus1 + 1;

sw = scat_matrix(out_angle_idx, inc_angle_idx);
ne = scat_matrix(out_angle_idx_plus1, inc_angle_idx_plus1);
se = scat_matrix(out_angle_idx, inc_angle_idx_plus1);
nw = scat_matrix(out_angle_idx_plus1, inc_angle_idx);

f1 = sw + (se - sw) .* inc_angle_frac;
f2 = nw + (ne - nw) .* inc_angle_frac;

val = f1 + (f2 - f1) .* out_angle_frac;

end
function alpha = fn_angle_from_probe_normal(ray_geometry, probe_coords)
% Simply computes the outgoing ray angle with respect to the probe normal.
% Will only be required if the probe angle is non-zero.
%
% INPUTS:
% - ray_geometry : array (2, 3)
%       3D coordinates of the start and end points of the first leg of the
%       ray (i.e. ray_geometry(1, :) is the 3D coords of the probe element
%       from which the ray comes from, ray_geometry(2, :) is the
%       coordinates of the first wall the ray interacts with).
% - probe_coords : array (probe_els, 3)
%       3D coordinates of the probe elements.
%
% OUTPUTS:
% - alpha : double
%       Outgoing angle the ray makes with the probe normal in radians.

norm_vector = zeros(3, 1);
norm_radius = 0;
ray_vector = zeros(3, 1);
ray_radius = 0;
dot_out = 0;

rot_mat = [-1, -1, 1];

for ii = 1:3
    norm_vector(ii, 1) = rot_mat(ii) * (probe_coords(end, 4-ii) - probe_coords(1, 4-ii));
    ray_vector(ii, 1) = ray_geometry(end, ii) - ray_geometry(1, ii);
    
    norm_radius = norm_radius + norm_vector(ii, 1) * norm_vector(ii, 1);
    ray_radius = ray_radius + ray_vector(ii, 1) * ray_vector(ii, 1);
    dot_out = dot_out + norm_vector(ii, 1) * ray_vector(ii, 1);
end

norm_radius = sqrt(norm_radius);
ray_radius = sqrt(ray_radius);

alpha = acos(dot_out/(norm_radius * ray_radius));

end
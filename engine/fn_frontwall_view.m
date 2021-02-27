function view = fn_frontwall_view(probe_coords, probe_angle, front_wall, mat_speeds, density, probe_frequency, probe_pitch)
% Computes the view arising from reflections from the front wall of the
% geometry. Note that this is an optional function, as these are not
% required in the contact case, or if the geometry is not being simulated.
%
% INPUTS:
% - probe_coords : array (probe_els, 3)
%       3D coordinates of the each probe element in the array.
% - probe_angle : double
%       Angle which the probe makes with respect to the geometry in
%       radians. Note that angle is zero when the probe is parallel with
%       the front wall of the geometry (and must always be zero in the
%       contact case).
% - front_wall : array (front_wall_points, 3)
%       3D coordinates of the front wall. Discretised into a number of
%       points.
% - mat_speeds : array (3, 1)
%       Speeds of the ray in each medium and for each mode. The order must
%       be [couplant_speed, solid_long_speed, solid_shear_speed]. Note that
%       all three must be provided, even in the contact case.
% - densities : array (2, 1)
%       Density of each medium which the ray passes through. The order must
%       be [couplant_density, solid_density]. Note that couplant density
%       must be provided, even in the contact case.
% - probe_freq : double
%       The frequency of the probe. Note that this is just the probe
%       frequency, not the frequency array which would be computed in the
%       multi-frequency model.
% - probe_pitch : double
%       The pitch of the probe array (i.e. the length of each probe element
%       in the plane being inspected).
%
% OUTPUTS:
% - View : struct (1, 1)
%       Structure containing the view generated from the frontwall. Only
%       one is returned as in the couplant, the only possible wave mode is
%       lontitudinal.

[probe_els, ~] = size(probe_coords);
[fw_points, ~] = size(front_wall);

probe_as_scatterer.image_block = probe_coords;
probe_as_scatterer.x_shape = 1;
probe_as_scatterer.z_shape = 1;
probe_as_scatterer.type = "image";

geometry = zeros(1, size(front_wall, 1), size(front_wall, 2));
geometry(1, :, :) = front_wall;

frontwall_path = fn_path_info( ...
    'Frontwall L-L', [1, 1], geometry, [mat_speeds(1), mat_speeds(1)], mat_speeds, [1], ...
    [0, 0], density, probe_frequency, probe_pitch, probe_coords ...
);

ray = fn_compute_ray(probe_as_scatterer, frontwall_path, 0);
view.name = 'Frontwall L-L';
view.min_times = zeros(probe_els ^ 2, 1);
view.probe_txrx = zeros(probe_els ^ 2, 2);
view.min_theta1 = zeros(probe_els ^ 2, 1);
view.min_theta2 = zeros(probe_els ^ 2, 1);
view.min_dist1 = zeros(probe_els ^ 2, 1);
view.bs1 = zeros(probe_els ^ 2, 1);
view.path_1 = ray;
view.tr1 = zeros(probe_els ^ 2, 1);
view.tau = ray.min_times;
view.directivity1 = zeros(probe_els ^ 2, 1);
view.directivity2 = zeros(probe_els ^ 2, 1);

el = 1;
for t_el = 1 : probe_els
    for r_el = 1 : probe_els
        
        view.min_times(el, 1) = ray.min_times(t_el, r_el);
        view.probe_txrx(el, :) = [t_el, r_el];
        view.min_theta1(el, 1) = ray.min_theta(t_el, r_el);
        view.min_theta2(el, 1) = ray.min_theta(r_el, t_el);
        view.min_dist1(el, 1) = ray.min_dists(t_el, r_el);
        view.bs1(el, 1) = ray.beam_spread(t_el, r_el);
        view.tr1(el, 1) = ray.tr_coeff(t_el, r_el);
        view.directivity1(el, 1) = ray.directivity(t_el, r_el);
        view.directivity2(el, 1) = ray.inv_directivity(t_el, r_el);
        
        el = el+1;
    end
end

end
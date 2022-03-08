function ray_weights = fn_compute_ray_weights(ray, freq_array)
% Computes ray weights from the ray path.
% INPUTS:
% - ray : output from fn_compute_ray function. Contains rays for which
%       weights will be computed.
% - freq_array : array
%       Non-zero frequencies, as output from fn_create_input_signal
%       function.
%
% OUTPUTS:
% - ray_weights : struct (1, 1)
%       Contains beamspread, directivity and transmission/reflection
%       coefficients along the forward and backward path for the ray.

% Relevant dimensions.
[probe_els, num_scatterers] = size(ray.min_times);
[num_npws, ~] = size(freq_array);

% Unpack path and scatterer info.
path_info = ray.path_info;
path_geometry = path_info.path_geometry;
if ~isstruct(path_geometry)
    no_walls = 0;
else
    [no_walls, ~, ~] = size(path_geometry);
end
speeds = path_info.speeds;
mat_speeds = path_info.material_speeds;
modes = path_info.modes;
probe_pitch = path_info.probe_pitch;
el_length = path_info.el_length;
scat_info = ray.scat_info;
probe_coords = path_info.probe_coords;
medium_ids = path_info.medium_ids;
density = path_info.densities;
scatterers = scat_info.image_block;
npw = path_info.npw;
[~, num_npws] = size(npw);
c44 = path_info.modulus / (2. * (1. + path_info.poisson));

% Initialise beamspread, directivity and transmission/reflection coeffs.
ray_weights.beamspread = zeros(probe_els, num_scatterers, num_npws);
ray_weights.inv_beamspread = zeros(probe_els, num_scatterers, num_npws);
ray_weights.directivity = zeros(probe_els, num_scatterers, num_npws);
ray_weights.transrefl = zeros(probe_els, num_scatterers, num_npws);
ray_weights.inv_transrefl = zeros(probe_els, num_scatterers, num_npws);
% Total weights for each frequency.
ray_weights.weights = zeros(probe_els, num_scatterers, num_npws);
ray_weights.inv_weights = zeros(probe_els, num_scatterers, num_npws);
% Initialise angles.
ray_weights.inc_theta = zeros(probe_els, num_scatterers, no_walls+1, 2);
ray_weights.out_theta = zeros(probe_els, num_scatterers, no_walls+1, 2);
ray_weights.inv_inc_theta = zeros(probe_els, num_scatterers, no_walls+1, 2);
ray_weights.inv_out_theta = zeros(probe_els, num_scatterers, no_walls+1, 2);
ray_weights.c_out = zeros(probe_els,1);
ray_weights.min_dists = zeros(probe_els, num_scatterers, no_walls+1, 4);

% If we are in the direct contact case
if ~isstruct(path_geometry)
    for scat = 1 : num_scatterers
        for tx = 1 : probe_els
            
            single_ray_leg_coords = zeros(2, 3);
            single_ray_leg_coords(1, :) = probe_coords(tx, :);
            single_ray_leg_coords(2, :) = scatterers(scat, :);
            inv_single_ray_geometry = flip(single_ray_leg_coords, 1);

            min_dists = fn_min_dists(single_ray_leg_coords);
            inv_min_dists = fn_min_dists(inv_single_ray_geometry);
                
            inc_out_angles = fn_inc_out_angles(min_dists, path_geometry);
            inc_out_angles(1, 2) = fn_angle_from_probe_normal(single_ray_leg_coords(1:2, :), probe_coords);
            inv_inc_out_angles = flip(inc_out_angles, 1);
            
            ray_weights.inc_theta(tx, scat, :, 1) = inc_out_angles(:, 1);
            ray_weights.inc_theta(tx, scat, :, 2) = inc_out_angles(:, 3);
            ray_weights.out_theta(tx, scat, :, 1) = inc_out_angles(:, 2);
            ray_weights.out_theta(tx, scat, :, 2) = inc_out_angles(:, 4);
            ray_weights.inv_inc_theta(tx, scat, :, 1) = inv_inc_out_angles(:, 1);
            ray_weights.inv_inc_theta(tx, scat, :, 2) = inv_inc_out_angles(:, 3);
            ray_weights.inv_out_theta(tx, scat, :, 1) = inv_inc_out_angles(:, 2);
            ray_weights.inv_out_theta(tx, scat, :, 2) = inv_inc_out_angles(:, 4);
            if scat==118
                ray_weights.c_out(tx) = inc_out_angles(end, 3);
            end
            ray_weights.min_dists(tx, scat, :, :) = min_dists;
            
            bs = fn_beamspread_2d(min_dists(:, 4), 0, speeds);
            ibs = fn_beamspread_2d(inv_min_dists(:, 4), 0, flip(speeds));
            
            tr = 1;
            
            dir = fn_sinc_directivity(inc_out_angles(1, 2), el_length, speeds(1)/freq_array(1));
            
            for npw_idx = 1 : num_npws
                
                % Note angles are not passed to beamspread functions:
                % because we are in direct contact case, there is only one
                % leg in the ray, and thus the virtual distance which is
                % calculated in this function is equal to the real
                % distance - no angles are required in the calculation.
                ray_weights.beamspread(tx, scat, npw_idx) = bs;
                ray_weights.inv_beamspread(tx, scat, npw_idx) = ibs;
                
                % As there are no transmissions or reflections, set the
                % coefficient equal to one.
                ray_weights.transrefl(tx, scat, npw_idx) = tr;
                ray_weights.inv_transrefl(tx, scat, npw_idx) = tr;
                
                % We are in direct contact path.
                % We want the outgoing angle from the probe.
                ray_weights.directivity(tx, scat, npw_idx) = dir;% * ...
%                     fn_linedir_lookup(inc_out_angles(1, 2), modes(1), min_dists(1, 4));
%                     fn_line_directivity(inc_out_angles(1, 2), mat_speeds(2)/freq_array(freq_idx), mat_speeds(3)/freq_array(freq_idx), modes(1), c44);
            end
        end
    end
    
% If we are in the skip contact case or immersion case, i.e. we
% have more than one leg.
else
    for scat = 1 : num_scatterers
        for tx = 1 : probe_els
            
            single_ray_leg_coords = zeros(no_walls+2, 3);
            single_ray_leg_coords(1, :) = probe_coords(tx, :);
            for wall = 1:no_walls
                single_ray_leg_coords(wall+1, :) = path_geometry(wall).coords(ray.wall_idxs(tx, scat, wall), :);
            end
            single_ray_leg_coords(no_walls+2, :) = scatterers(scat, :);
            inv_single_ray_geometry = flip(single_ray_leg_coords, 1);

            % Get the minimum distances.
            min_dists = fn_min_dists(single_ray_leg_coords);
            inv_min_dists = fn_min_dists(inv_single_ray_geometry);

            % Get the angles.
            inc_out_angles = fn_inc_out_angles(min_dists, path_geometry);
            inc_out_angles(1, 2) = fn_angle_from_probe_normal(single_ray_leg_coords(1:2, :), probe_coords);
            inv_inc_out_angles = flip(inc_out_angles, 1);
            
            ray_weights.inc_theta(tx, scat, :, 1) = inc_out_angles(:, 1);
            ray_weights.inc_theta(tx, scat, :, 2) = inc_out_angles(:, 3);
            ray_weights.out_theta(tx, scat, :, 1) = inc_out_angles(:, 2);
            ray_weights.out_theta(tx, scat, :, 2) = inc_out_angles(:, 4);
            ray_weights.inv_inc_theta(tx, scat, :, 1) = inv_inc_out_angles(:, 1);
            ray_weights.inv_inc_theta(tx, scat, :, 2) = inv_inc_out_angles(:, 3);
            ray_weights.inv_out_theta(tx, scat, :, 1) = inv_inc_out_angles(:, 2);
            ray_weights.inv_out_theta(tx, scat, :, 2) = inv_inc_out_angles(:, 4);
            if scat==118
                ray_weights.c_out(tx) = inc_out_angles(end, 3);
            end
            ray_weights.min_dists(tx, scat, :, :) = min_dists;
            
            bs = fn_beamspread_2d(min_dists(:, 4), inc_out_angles(:, 1), speeds);
            ibs = fn_beamspread_2d(inv_min_dists(:, 4), inv_inc_out_angles(:, 2), flip(speeds));
            
            inv_angles = conj(asin(sin(inc_out_angles(1:end-1, 1)) .* speeds(2:end) ./ speeds(1:end-1)));
            tr = fn_TR_coeff(medium_ids, modes, inc_out_angles(1:end-1, 1), mat_speeds, density);
            itr = fn_TR_coeff(flip(medium_ids), flip(modes), inv_angles, mat_speeds, density);
            
            if medium_ids(1) == 1
                dir = fn_sinc_directivity(inc_out_angles(1, 2), el_length, speeds(1)/freq_array(1));
            else
                dir = fn_sinc_directivity(inc_out_angles(1, 2), el_length, speeds(end)/freq_array(1));
            end
            
            for npw_idx = 1 : num_npws
                
                % Beamspread
                ray_weights.beamspread(tx, scat, npw_idx) = bs;
                ray_weights.inv_beamspread(tx, scat, npw_idx) = ibs;
                
                % Trans/refl
                ray_weights.transrefl(tx, scat, npw_idx) = tr;
                ray_weights.inv_transrefl(tx, scat, npw_idx) = itr;

                % We are in skip contact or immersion case. Check which,
                % and then compute directivity.
                ray_weights.directivity(tx, scat, npw_idx) = dir;
            end
        end
    end
end

theta_for_dir = reshape(ray_weights.out_theta(:, :, 1, 1), [], 1);
modes_for_dir = modes(1);
dists_for_dir = reshape(ray_weights.min_dists(:, :, 1, 4), [], 1);

useabs = 0;
for npw_idx = 1:size(npw, 2)
    if useabs
        linedir = abs(fn_linedir_lookup(theta_for_dir, modes_for_dir, dists_for_dir, npw(npw_idx)));
    else
        linedir = fn_linedir_lookup(theta_for_dir, modes_for_dir, dists_for_dir, npw(npw_idx));
    end
    ray_weights.directivity(:, :, npw_idx) = ray_weights.directivity(:, :, npw_idx) .* reshape(linedir, probe_els, num_scatterers); %
end

% Calculate the total ray weights, for convenience in calculation later.
ray_weights.weights = ( ...
    ray_weights.transrefl .* ray_weights.directivity .* ray_weights.beamspread ...
);
ray_weights.inv_weights = ( ...
    ray_weights.inv_transrefl .* ray_weights.inv_beamspread .* ray_weights.directivity ...
) .* sqrt(speeds(end)/reshape(freq_array', 1,1,1));
end
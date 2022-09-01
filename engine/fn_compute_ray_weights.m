function ray_weights = fn_compute_ray_weights(ray, freq_array, varargin)
% Computes ray weights from the ray path.
% INPUTS:
% - ray : output from fn_compute_ray function. Contains rays for which
%       weights will be computed.
% - freq_array : array
%       Non-zero frequencies, as output from fn_create_input_signal
%       function.
% - all_geometry : struct : OPTIONAL
%       Optionally pass in to make plotting the ray possible.
%
% OUTPUTS:
% - ray_weights : struct (1, 1)
%       Contains beamspread, directivity and transmission/reflection
%       coefficients along the forward and backward path for the ray.

debug = false;

% Relevant dimensions.
[probe_els, num_scatterers] = size(ray.min_times);
[num_freqs, ~] = size(freq_array);

for arg = 1:nargin-2
    argument = varargin{arg};
    if isa(argument, 'struct')
        all_geometry = argument;
    end
end

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
scatterers = [scat_info.x, scat_info.y, scat_info.z];
npw = path_info.npw;
c44 = path_info.modulus / (2. * (1. + path_info.poisson));

% Initialise beamspread, directivity and transmission/reflection coeffs.
ray_weights.beamspread = zeros(probe_els, num_scatterers, num_freqs);
ray_weights.inv_beamspread = zeros(probe_els, num_scatterers, num_freqs);
ray_weights.directivity = zeros(probe_els, num_scatterers, num_freqs);
ray_weights.idirectivity = zeros(probe_els, num_scatterers, num_freqs);
ray_weights.transrefl = zeros(probe_els, num_scatterers, num_freqs);
ray_weights.inv_transrefl = zeros(probe_els, num_scatterers, num_freqs);
% Total weights for each frequency.
ray_weights.weights = zeros(probe_els, num_scatterers, num_freqs);
ray_weights.inv_weights = zeros(probe_els, num_scatterers, num_freqs);
% Initialise angles.
% ray_weights.inc_theta = zeros(probe_els, num_scatterers, no_walls+1, 2);
% ray_weights.out_theta = zeros(probe_els, num_scatterers, no_walls+1, 2);
% ray_weights.inv_inc_theta = zeros(probe_els, num_scatterers, no_walls+1, 2);
% ray_weights.inv_out_theta = zeros(probe_els, num_scatterers, no_walls+1, 2);
% ray_weights.c_out = zeros(probe_els,1);
ray_weights.min_dists = zeros(probe_els, num_scatterers, no_walls+1, 4);

% if debug
%     close all
%     figure(3)
%     hold on
%     for ii = 1:size(all_geometry, 1)
%         plot([all_geometry(ii).point1(1), all_geometry(ii).point2(1)], -[all_geometry(ii).point1(3), all_geometry(ii).point2(3)], 'r')
%     end
% end

% If we are in the direct contact case
if ~isstruct(path_geometry)
    for tx = 1 : probe_els

        single_ray_leg_coords = permute(squeeze(ray.coords(tx, :, :, :)), [2,3,1]);
        inv_single_ray_geometry = flip(single_ray_leg_coords, 1);

        min_dists = fn_min_dists(single_ray_leg_coords);
        inv_min_dists = fn_min_dists(inv_single_ray_geometry);

        inc_out_angles = fn_inc_out_angles(min_dists, path_geometry);
        inc_out_angles(1, 2) = fn_angle_from_probe_normal(single_ray_leg_coords(1:2, :), probe_coords);
        inv_inc_out_angles = flip(inc_out_angles, 1);            
%         inv_inc_out_angles = fn_inc_out_angles(inv_min_dists, path_geometry);
% 
%         ray_weights.inc_theta(tx, scat, :, 1) = inc_out_angles(:, 1);
%         ray_weights.inc_theta(tx, scat, :, 2) = inc_out_angles(:, 3);
%         ray_weights.out_theta(tx, scat, :, 1) = inc_out_angles(:, 2);
%         ray_weights.out_theta(tx, scat, :, 2) = inc_out_angles(:, 4);
%         ray_weights.inv_inc_theta(tx, scat, :, 1) = inv_inc_out_angles(:, 1);
%         ray_weights.inv_inc_theta(tx, scat, :, 2) = inv_inc_out_angles(:, 3);
%         ray_weights.inv_out_theta(tx, scat, :, 1) = inv_inc_out_angles(:, 2);
%         ray_weights.inv_out_theta(tx, scat, :, 2) = inv_inc_out_angles(:, 4);
%         if scat==118
%             ray_weights.c_out(tx) = inc_out_angles(end, 3);
%         end
        ray_weights.min_dists(tx, :, :, :) = permute(min_dists, [3, 1, 2]);

%         cos_sin = [cos(inc_out_angles(1, 2)), sin(inc_out_angles(1, 2))];
            
%         if debug
%             plot(single_ray_leg_coords(:, 1), -single_ray_leg_coords(:, 3), 'Color', [.5,.5,.5])
%         end
        for freq_idx = 1 : num_freqs
            % Note angles are not passed to beamspread functions:
            % because we are in direct contact case, there is only one
            % leg in the ray, and thus the virtual distance which is
            % calculated in this function is equal to the real
            % distance - no angles are required in the calculation.
            ray_weights.beamspread(tx, :, freq_idx) = fn_beamspread_2d(reshape(min_dists(:, 4, :), no_walls+1, num_scatterers), 0, speeds, freq_array(freq_idx));
            ray_weights.inv_beamspread(tx, :, freq_idx) = fn_beamspread_2d(reshape(inv_min_dists(:, 4, :), no_walls+1, num_scatterers), 0, flip(speeds), freq_array(freq_idx));

            % As there are no transmissions or reflections, set the
            % coefficient equal to one.
            ray_weights.transrefl(tx, :, freq_idx) = 1;
            ray_weights.inv_transrefl(tx, :, freq_idx) = 1;
            
            % We are in direct contact path.
            % We want the outgoing angle from the probe.
            ray_weights.directivity(tx, :, freq_idx) = ...
                fn_sinc_directivity(reshape(inc_out_angles(1, 2, :), 1, num_scatterers), el_length, speeds(1)/freq_array(freq_idx));
            ray_weights.idirectivity(tx, :, freq_idx) = ...
                fn_sinc_directivity(reshape(inc_out_angles(1, 2, :), 1, num_scatterers), el_length, speeds(1)/freq_array(freq_idx));% * ...
%                     cos_sin(modes(1)+1);
%                     fn_line_directivity(inc_out_angles(1, 2), mat_speeds(2)/freq_array(freq_idx), mat_speeds(3)/freq_array(freq_idx), modes(1));
            if npw == 0
                ray_weights.directivity(tx, :, freq_idx) = ray_weights.directivity(tx, :, freq_idx) .* ...
                    fn_line_directivity(reshape(inc_out_angles(1, 2, :), 1, num_scatterers), mat_speeds(2)/freq_array(freq_idx), mat_speeds(3)/freq_array(freq_idx), modes(1));
            end
        end


%         for scat = 1 : num_scatterers
%             for freq_idx = 1 : num_freqs
%                 
%                 % Note angles are not passed to beamspread functions:
%                 % because we are in direct contact case, there is only one
%                 % leg in the ray, and thus the virtual distance which is
%                 % calculated in this function is equal to the real
%                 % distance - no angles are required in the calculation.
%                 ray_weights.beamspread(tx, scat, freq_idx) = fn_beamspread_2d(min_dists(:, 4, scat), 0, speeds, freq_array(freq_idx));
%                 ray_weights.inv_beamspread(tx, scat, freq_idx) = fn_beamspread_2d(inv_min_dists(:, 4, scat), 0, flip(speeds), freq_array(freq_idx));
%                 
%                 % As there are no transmissions or reflections, set the
%                 % coefficient equal to one.
%                 ray_weights.transrefl(tx, scat, freq_idx) = 1;
%                 ray_weights.inv_transrefl(tx, scat, freq_idx) = 1;
%                 
%                 % We are in direct contact path.
%                 % We want the outgoing angle from the probe.
%                 ray_weights.directivity(tx, scat, freq_idx) = ...
%                     fn_sinc_directivity(inc_out_angles(1, 2, scat), el_length, speeds(1)/freq_array(freq_idx));
%                 ray_weights.idirectivity(tx, scat, freq_idx) = ...
%                     fn_sinc_directivity(inc_out_angles(1, 2, scat), el_length, speeds(1)/freq_array(freq_idx));% * ...
% %                     cos_sin(modes(1)+1);
% %                     fn_line_directivity(inc_out_angles(1, 2), mat_speeds(2)/freq_array(freq_idx), mat_speeds(3)/freq_array(freq_idx), modes(1));
%                 if npw == 0
%                     ray_weights.directivity(tx, scat, freq_idx) = ray_weights.directivity(tx, scat, freq_idx) * ...
%                         fn_line_directivity(inc_out_angles(1, 2, scat), mat_speeds(2)/freq_array(freq_idx), mat_speeds(3)/freq_array(freq_idx), modes(1));
%                 end
%             end
%         end
    end
    
% If we are in the skip contact case or immersion case, i.e. we
% have more than one leg.
else
    for tx = 1 : probe_els
        
        single_ray_leg_coords = permute(squeeze(ray.coords(tx, :, :, :)), [2,3,1]);
        inv_single_ray_geometry = flip(single_ray_leg_coords, 1);

        min_dists = fn_min_dists(single_ray_leg_coords);
        inv_min_dists = fn_min_dists(inv_single_ray_geometry);

        inc_out_angles = fn_inc_out_angles(min_dists, path_geometry);
        inc_out_angles(1, 2) = fn_angle_from_probe_normal(single_ray_leg_coords(1:2, :), probe_coords);
        inv_inc_out_angles = flip(inc_out_angles, 1);   
%         inv_inc_out_angles = fn_inc_out_angles(inv_min_dists, path_geometry);
% 
%         ray_weights.inc_theta(tx, scat, :, 1) = inc_out_angles(:, 1);
%         ray_weights.inc_theta(tx, scat, :, 2) = inc_out_angles(:, 3);
%         ray_weights.out_theta(tx, scat, :, 1) = inc_out_angles(:, 2);
%         ray_weights.out_theta(tx, scat, :, 2) = inc_out_angles(:, 4);
%         ray_weights.inv_inc_theta(tx, scat, :, 1) = inv_inc_out_angles(:, 1);
%         ray_weights.inv_inc_theta(tx, scat, :, 2) = inv_inc_out_angles(:, 3);
%         ray_weights.inv_out_theta(tx, scat, :, 1) = inv_inc_out_angles(:, 2);
%         ray_weights.inv_out_theta(tx, scat, :, 2) = inv_inc_out_angles(:, 4);
%         if scat==118
%             ray_weights.c_out(tx) = inc_out_angles(end, 3);
%         end
        ray_weights.min_dists(tx, :, :, :) = permute(min_dists, [3, 1, 2]);

%         cos_sin = [cos(inc_out_angles(1, 2)), sin(inc_out_angles(1, 2))];
            
%         if debug
%             plot(single_ray_leg_coords(:, 1), -single_ray_leg_coords(:, 3), 'Color', [.5,.5,.5])
%         end
        for freq_idx = 1 : num_freqs
            % Note angles are not passed to beamspread functions:
            % because we are in direct contact case, there is only one
            % leg in the ray, and thus the virtual distance which is
            % calculated in this function is equal to the real
            % distance - no angles are required in the calculation.            
            ray_weights.beamspread(tx, :, freq_idx) = fn_beamspread_2d(reshape(min_dists(:, 4, :), no_walls+1, num_scatterers), reshape(inc_out_angles(:, 1, :), no_walls+1, num_scatterers), speeds, freq_array(freq_idx));
            ray_weights.inv_beamspread(tx, :, freq_idx) = fn_beamspread_2d(reshape(inv_min_dists(:, 4, :), no_walls+1, num_scatterers), reshape(inv_inc_out_angles(:, 2, :), no_walls+1, num_scatterers), flip(speeds), freq_array(freq_idx));

            % Trans/refl      
            inv_angles = reshape(conj(asin(sin(inc_out_angles(1:end-1, 1, :)) .* speeds(2:end) ./ speeds(1:end-1))), no_walls, num_scatterers);
            ray_weights.transrefl(tx, :, freq_idx) = fn_TR_coeff(medium_ids, modes, reshape(inc_out_angles(1:end-1, 1, :), no_walls, num_scatterers), mat_speeds, density);
            ray_weights.inv_transrefl(tx, :, freq_idx) = fn_TR_coeff(flip(medium_ids), flip(modes), inv_angles, mat_speeds, density);

            % We are in skip contact or immersion case. Check which,
            % and then compute directivity.
            if medium_ids(1) == 1 % If solid then contact.
                ray_weights.directivity(tx, :, freq_idx) = ...
                    fn_sinc_directivity(inc_out_angles(1, 2, :), el_length, speeds(1)/freq_array(freq_idx));
                ray_weights.idirectivity(tx, :, freq_idx) = ...
                    fn_sinc_directivity(inc_out_angles(1, 2, :), el_length, speeds(1)/freq_array(freq_idx));% * ...
%                         cos_sin(modes(1)+1);
%                         fn_line_directivity(inc_out_angles(1, 2), mat_speeds(2)/freq_array(freq_idx), mat_speeds(3)/freq_array(freq_idx), modes(1));
                if npw == 0
                    ray_weights.directivity(tx, :, freq_idx) = ray_weights.directivity(tx, :, freq_idx) .* ...
                        fn_line_directivity(reshape(inc_out_angles(1, 2, :), 1, num_scatterers), mat_speeds(2)/freq_array(freq_idx), mat_speeds(3)/freq_array(freq_idx), modes(1));
                end
                
            else % Must be liquid, so immersion.
                ray_weights.directivity(tx, :, freq_idx) = fn_sinc_directivity(inc_out_angles(1, 2, :), el_length, speeds(end)/freq_array(freq_idx));
            end
        end
            
%         for scat = 1 : num_scatterers
%             for freq_idx = 1 : num_freqs
%                 
%                 % Beamspread
%                 ray_weights.beamspread(tx, scat, freq_idx) = fn_beamspread_2d(min_dists(:, 4, scat), inc_out_angles(:, 1, scat), speeds, freq_array(freq_idx));
%                 ray_weights.inv_beamspread(tx, scat, freq_idx) = fn_beamspread_2d(inv_min_dists(:, 4, scat), inv_inc_out_angles(:, 2, scat), flip(speeds), freq_array(freq_idx));
%                 
%                 % Trans/refl
%                 inv_angles = conj(asin(sin(inc_out_angles(1:end-1, 1, scat)) .* speeds(2:end) ./ speeds(1:end-1)));
%                 ray_weights.transrefl(tx, scat, freq_idx) = fn_TR_coeff(medium_ids, modes, inc_out_angles(1:end-1, 1, scat), mat_speeds, density);
%                 ray_weights.inv_transrefl(tx, scat, freq_idx) = fn_TR_coeff(flip(medium_ids), flip(modes), inv_angles, mat_speeds, density);
% 
%                 % We are in skip contact or immersion case. Check which,
%                 % and then compute directivity.
%                 if medium_ids(1) == 1 % If solid then contact.
%                     ray_weights.directivity(tx, scat, freq_idx) = ...
%                         fn_sinc_directivity(inc_out_angles(1, 2, scat), el_length, speeds(1)/freq_array(freq_idx));
%                     ray_weights.idirectivity(tx, scat, freq_idx) = ...
%                         fn_sinc_directivity(inc_out_angles(1, 2, scat), el_length, speeds(1)/freq_array(freq_idx));% * ...
% %                         cos_sin(modes(1)+1);
% %                         fn_line_directivity(inc_out_angles(1, 2), mat_speeds(2)/freq_array(freq_idx), mat_speeds(3)/freq_array(freq_idx), modes(1));
%                     if npw == 0
%                         ray_weights.directivity(tx, scat, freq_idx) = ray_weights.directivity(tx, scat, freq_idx) * ...
%                             fn_line_directivity(inc_out_angles(1, 2, scat), mat_speeds(2)/freq_array(freq_idx), mat_speeds(3)/freq_array(freq_idx), modes(1));
%                     end
%                 else % Must be liquid, so immersion.
%                     ray_weights.directivity(tx, scat, freq_idx) = fn_sinc_directivity(inc_out_angles(1, 2, scat), el_length, speeds(end)/freq_array(freq_idx));
%                 end
%             end
%         end
    end
end

if npw ~= 0
    theta_for_dir = reshape(ray_weights.out_theta(:, :, 1, 1), [], 1);
    theta_for_dir = mod(pi/2 + reshape(ray_weights.out_theta(:, :, 1, 1), [], 1), pi) - pi/2;
    modes_for_dir = modes(1);
    dists_for_dir = reshape(ray_weights.min_dists(:, :, 1, 4), [], 1);

    useabs = 0;
    linedir = fn_linedir_lookup(theta_for_dir, modes_for_dir, dists_for_dir, npw);
    if useabs
        linedir = abs(linedir);
    end
    ray_weights.directivity(:, :, 1) = ray_weights.directivity(:, :, 1) .* reshape(linedir, probe_els, num_scatterers);
end


% Calculate the total ray weights, for convenience in calculation later.
% Note sqrt(lambda) here comes from the fact that the scattering matrix is
% defined as the amplitude of a scattered wave that would be measured if
% the distance between defect and receiver is one wavelength (Zhang, 2008)
ray_weights.weights = ( ...
    ray_weights.transrefl .* ray_weights.directivity .* ray_weights.beamspread ...
) ./ c44;% ./ (density(medium_ids(1)+1) * speeds(1));
ray_weights.inv_weights = ( ...
    ray_weights.inv_transrefl .* ray_weights.inv_beamspread .* ray_weights.directivity ...
) .* sqrt(speeds(end)/reshape(freq_array', 1,1,num_freqs)) .* (1i * density(medium_ids(1)+1) * speeds(1));

% Turn off weight breakdown - not necessary to keep them unless debugging.
if ~debug
    ray_weights = rmfield(ray_weights, 'min_dists');
    ray_weights = rmfield(ray_weights, 'beamspread');
    ray_weights = rmfield(ray_weights, 'inv_beamspread');
    ray_weights = rmfield(ray_weights, 'directivity');
    ray_weights = rmfield(ray_weights, 'idirectivity');
    ray_weights = rmfield(ray_weights, 'transrefl');
    ray_weights = rmfield(ray_weights, 'inv_transrefl');
end

end
function Views = fn_geometry_views(probe_coords, all_geometries, mat_speeds, densities, probe_freq, probe_pitch, max_no_refl)
% Computes signals which come from reflections of the wall geometry
% (including frontwall reflections).
%
% INPUTS:
% - probe_coords : array (probe_els, 3)
%       3D coordinates of the each probe element in the array.
% - mat_speeds : array (3, 1)
%       Speeds of the ray in each medium and for each mode. The order must
%       be [couplant_speed, solid_long_speed, solid_shear_speed]. Note that
%       all three must be provided, even in the contact case.
%
% OUTPUTS:
% - Views : struct (num_geom_views, 1)
%       Structure containing all of the views arising from the geometry
%       walls. For a rectangular block in immersion, there will be 5
%       (frontwall, backwall L-L, bw L-T, bw T-L, bw T-T).

[probe_els, ~] = size(probe_coords);

% Set up variables required for calculations.
probe_as_scatterer.image_block = probe_coords;
probe_as_scatterer.x_shape = 1;
probe_as_scatterer.z_shape = 1;
probe_as_scatterer.type = "image";

couplant_spd = mat_speeds(1);
mat_long_spd = mat_speeds(2);
mat_shear_spd = mat_speeds(3);

num_walls = size(all_geometries, 1);

wall_names = repmat("", num_walls, 1);
for wall = 1:size(all_geometries, 1)
    wall_names(wall) = all_geometries(wall).name;
end

view.name = "name";
view.min_times = zeros(probe_els ^ 2, 1);
view.probe_txrx = zeros(probe_els ^ 2, 2);
view.valid_paths = zeros(probe_els ^ 2, 1);
view.weights = zeros(probe_els^2, 1);

mat_1_spd = [mat_long_spd, mat_long_spd, mat_shear_spd, mat_shear_spd];
mat_2_spd = [mat_long_spd, mat_shear_spd, mat_long_spd, mat_shear_spd];
modes = [[0, 0]; [0, 1]; [1, 0]; [1, 1]];
names = [["L", "L"]; ["L", "T"]; ["T", "L"]; ["T", "T"]];

% If we are in contact.
if ~ismember("F", wall_names)
    Views = repmat(view, 4*num_walls, 1);
    
    view_speeds = zeros(4, 2); % Shape (no_views, no_legs)
    view_modes  = zeros(4, 2);
    
    for path = 1:4
        view_speeds(path, 1) = mat_1_spd(path);
        view_speeds(path, 2) = mat_2_spd(path);
        view_modes(path, 1)  = modes(path, 1);
        view_modes(path, 2)  = modes(path, 2);
    end
    
    walls = [2];
    medium_ids = [1, 1];
    
    non_fw_geometries = all_geometries;
    
% If we are in immersion.
else
    Views = repmat(view, 4*(num_walls-1)+1, 1);
    num_walls = num_walls - 1;
    
    view_speeds = zeros(4, 4); % Shape (no_views, no_legs)
    view_modes  = zeros(4, 4);
    
    for view = 1:4
        view_speeds(view, :) = [couplant_spd, mat_1_spd(view), mat_2_spd(view), couplant_spd];
        view_modes(view, :)  = [1, modes(view, 1), modes(view, 2), 1];
    end
    
    walls = [1, 2, 1];
    medium_ids = [0, 1, 1, 0];
    
    where_F = logical(wall_names=="F");
    
    non_fw_geometries = repmat(all_geometries(where_F), num_walls, 2);
    non_fw_geometries(:, 2) = all_geometries(~where_F);
    
end
    
% Loop over all of the walls.
for wall = 1:num_walls
    % Loop over modes for single reflection walls.
    for mode = 1:4
        view_idx = 4*(wall-1) + mode;

        % Get the path info.
        path_info = fn_path_info( ...
            sprintf("%s %s %s", names(mode, 1), non_fw_geometries(wall).name, names(mode, 2)), ...
            sprintf("%s %s %s", names(mode, 2), non_fw_geometries(wall).name, names(mode, 1)), ...
            view_modes(mode, :), ...
            non_fw_geometries(wall, :), ...
            view_speeds(mode, :), ...
            mat_speeds, ...
            walls, ...
            medium_ids, ...
            densities, ...
            probe_freq, ...
            probe_pitch, ...
            probe_coords ...
        );

        % Compute the ray.
        ray = fn_compute_ray(probe_as_scatterer, path_info, all_geometries, probe_freq);

        % Initialise this view's arrays.
        Views(view_idx).name = path_info.name;
        Views(view_idx).min_times = zeros(probe_els ^ 2, 1);
        Views(view_idx).probe_txrx = zeros(probe_els ^ 2, 2);
        Views(view_idx).valid_paths = zeros(probe_els ^ 2, 1);
        Views(view_idx).ray = ray;

        el = 1;
        % Assemble view from ray.
        for t_el = 1 : probe_els
            for r_el = 1 : probe_els

                Views(view_idx).min_times(el, 1) = ray.min_times(t_el, r_el);
                Views(view_idx).probe_txrx(el, :) = [t_el, r_el];
                Views(view_idx).valid_paths(el) = ray.valid_paths(t_el) .* ray.valid_paths(r_el);

                el = el+1;
            end
        end

        % Get the ray weights.
        [~, ~, num_freqs] = size(ray.weights.weights);
        Views(view_idx).weights = zeros(probe_els^2, num_freqs);
        el = 1;
        for ii = 1 : probe_els
            for jj = 1 : probe_els
                Views(view_idx).weights(el, :) = ( ...
                    ray.weights.weights(ii, jj, :) * ray.weights.inv_directivity(jj, ii, :) ...
                );

            el = el+1;
            end
        end


    end
end
    
end
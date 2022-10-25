function Views = fn_make_geometry_views(probe_coords, all_geometries, mat_speeds, densities, probe_freq, probe_pitch, el_length, max_no_refl, npw, is_contact, is_frontwall)
% Computes signals which come from reflections of the wall geometry
% (including frontwall reflections).
%
% INPUTS:
% - probe_coords : array (probe_els, 3)
%       3D coordinates of the each probe element in the array.
% - all_geometries : struct (total_num_walls, 1)
%       Structures containing all geometries. The only ones of interest are
%       backwalls ('Bj' for the jth wall) in contact, and the back and
%       frontwalls ('F') in immersion.
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
% - max_no_refl : integer
%       The maximum number of reflections to model from various backwalls.
%       Note that in the immersion case, this does not include the two
%       frontwall legs (i.e. probe -(L)-> F -(X)-> B -(X)-> F -(L)-> probe
%       only counts as one reflection). At the moment only supports two
%       reflections.
%
% OUTPUTS:
% - Views : struct (num_geom_views, 1)
%       Structure containing all of the views arising from the geometry
%       walls. For a rectangular block in immersion, there will be 5
%       (frontwall, backwall L-L, bw L-T, bw T-L, bw T-T).

[probe_els, ~] = size(probe_coords);

% Set up variables required for calculations.
probe_as_scatterer.x = probe_coords(:, 1);
probe_as_scatterer.y = probe_coords(:, 2);
probe_as_scatterer.z = probe_coords(:, 3);
% probe_as_scatterer.x_shape = 1;
% probe_as_scatterer.z_shape = 1;
probe_as_scatterer.type = "image";

couplant_spd = mat_speeds(1);

num_walls = size(all_geometries, 1);

wall_names = repmat("", num_walls, 1);
for wall = 1:size(all_geometries, 1)
    wall_names(wall) = all_geometries(wall).name;
end

names = ["L", "T"];

if ~is_frontwall
%% If there isn't a frontwall
    tot_num_views = 0;
    for refl = 1:max_no_refl
        tot_num_views = tot_num_views + ...
            2^(refl + 1) * num_walls * (num_walls - 1)^(refl - 1);
    end
    
    view_el = 1;
    for refl = 1:max_no_refl
        
        % Generate mode indices
        char_modes = dec2bin([0:2^(refl+1)-1]);
        view_modes = zeros(size(char_modes));
        for col = 1:size(char_modes, 2)
            view_modes(:, col) = str2num(char_modes(:, col)); %#ok<*ST2NM>
        end
        view_speeds = mat_speeds(view_modes+2);

        medium_ids = ones(refl+1, 1);
        
        % Generate wall indices
        if num_walls ~= 1
            wall_idxs = dec2base([0:num_walls^refl-1], num_walls);
        else
            % If there is only one wall, then dec2base does not work, but
            % only one wall so prescribe it as 0.
            wall_idxs = '0';
        end
        wall_idxs1 = zeros(size(wall_idxs));
        for col = 1:size(wall_idxs, 2)
            wall_idxs1(:, col) = str2num(wall_idxs(:, col)) + 1;
        end
        
        % Get the valid ones.
        wall_idxs = zeros(size(wall_idxs));
        this_row = 1;
        for row = 1:size(wall_idxs, 1)
            is_adjacent = 0;
            for col = 1:size(wall_idxs, 2)-1
                if wall_idxs1(row, col) == wall_idxs1(row, col+1)
                    is_adjacent = 1;
                    break
                end
            end
            if ~is_adjacent
                wall_idxs(this_row, :) = wall_idxs1(row, :);
                this_row = this_row+1;
            end
        end
        wall_idxs(wall_idxs(:, 1) == 0, :) = [];

        view_geometries = all_geometries(wall_idxs);

        for wall = 1:size(view_geometries, 1)
            for view = 1:size(view_modes, 1)
                % Get the name of this path.
                name = sprintf("%s", names(view_modes(view, 1)+1));
                rev_name = sprintf("%s", names(view_modes(view, end)+1));
                for mode = 2:size(view_modes, 2)
                    name = sprintf("%s %s %s", name, view_geometries(wall, mode-1).name, names(view_modes(view, mode)+1));
                    rev_name = sprintf("%s %s %s", rev_name, view_geometries(wall, end+2-mode).name, names(view_modes(view, end+1-mode)+1));
                end
                
                path_info = fn_path_info( ...
                    name, ...
                    rev_name, ...
                    view_modes(view, :), ...
                    view_geometries(wall, :).', ...
                    view_speeds(view, :), ...
                    mat_speeds, ...
                    medium_ids, ...
                    densities, ...
                    probe_freq, ...
                    probe_pitch, ...
                    el_length, ...
                    probe_coords, ...
                    npw ...
                );
                
                if exist("Views", "var")==0
                    Views = repmat(fn_create_geometry_view(path_info, all_geometries, probe_as_scatterer), tot_num_views, 1);
                else
                    Views(view_el) = fn_create_geometry_view(path_info, all_geometries, probe_as_scatterer);
                end
                view_el = view_el + 1;
            end
        end
    end
    
else
%% If we there is a front wall
    
    where_F = logical(wall_names=="F");
    idx_F = find(where_F);
    
    for refl = 1:max_no_refl
        % Generate mode indices
        char_modes = dec2bin(0:2^(refl+1)-1);
        view_modes = zeros(size(char_modes));
        for col = 1:size(char_modes, 2)
            view_modes(:, col) = str2num(char_modes(:, col)); %#ok<*ST2NM>
        end
        
        wall_idxs1 = dec2base(0:num_walls^refl-1, num_walls);
        wall_idxs = zeros(size(wall_idxs1));
        for col = 1:size(wall_idxs1, 2)
            wall_idxs(:, col) = str2num(wall_idxs1(:, col)) + 1;
        end
        
        % If immersion, we'll add the front wall to the beginning and end
        % later. Rays cannot join a wall to itself (i.e. a wall cannot be
        % adjacent to itself within a path), so remove all paths which
        % start or end with the front wall.
        wall_idxs = wall_idxs(and(wall_idxs(:, 1) ~= idx_F, wall_idxs(:, end) ~= idx_F), :);
        
        % Check for self-adjacency of all walls within paths.
        for row = 1:size(wall_idxs, 2)-1
            wall_idxs = wall_idxs(wall_idxs(:, row) ~= wall_idxs(:, row+1), :);
        end
        
        num_views_in_refl = size(wall_idxs, 1) * size(view_modes, 1);
        
        % If contact, do not need to add any walls to start or end.
        if is_contact
            view_speeds = mat_speeds(view_modes+2);
            view_geometries = all_geometries(wall_idxs);
            medium_ids = ones(refl+1, 1);
        % If immersion, add front wall to start and end.
        else
            view_speeds = cat(2, mat_speeds(1)*ones(size(view_modes, 1), 1), mat_speeds(view_modes+2), mat_speeds(1)*ones(size(view_modes, 1), 1));
            wall_idxs = cat(2, idx_F*ones(size(wall_idxs, 1), 1), wall_idxs, idx_F*ones(size(wall_idxs, 1), 1));
            view_geometries = all_geometries(wall_idxs);
            medium_ids = ones(refl+3, 1);
            medium_ids(1) = 0;
            medium_ids(end) = 0;
        end
        
        view_el = 1;
        for wall = 1:size(view_geometries, 1)
            for view = 1:size(view_modes, 1)
                name = sprintf('%s', names(view_modes(view, 1)+1));
                rev_name = sprintf('%s', names(view_modes(view, end)+1));
                for leg = 2:size(view_modes, 2)
                    name = sprintf('%s %s %s', name, wall_names(wall_idxs(wall, leg-1)), names(view_modes(view, leg)+1));
                    rev_name = sprintf('%s %s %s', rev_name, wall_names(wall_idxs(wall, end+2-leg)), names(view_modes(view, end+1-leg)+1));
                end
                
                path_info = fn_path_info( ...
                    name, ...
                    rev_name, ...
                    view_modes(view, :), ...
                    view_geometries(:, wall), ...
                    view_speeds(view, :), ...
                    mat_speeds, ...
                    medium_ids, ...
                    densities, ...
                    probe_freq, ...
                    probe_pitch, ...
                    el_length, ...
                    probe_coords, ...
                    npw ...
                );
                
                if exist("refl_Views", "var")==0
                    refl_Views = repmat(fn_create_geometry_view(path_info, all_geometries, probe_as_scatterer), num_views_in_refl, 1);
                else
                    refl_Views(view_el) = fn_create_geometry_view(path_info, all_geometries, probe_as_scatterer);
                end
                view_el = view_el + 1;
            end
        end
        if exist("Views", "var")==0
            Views = refl_Views;
        % Difficult to work out how big Views should be ahead of time due
        % to self-adjacency of walls for generalised number of reflections.
        % As such, precompute for each reflection and then expand Views
        % each time.
        else
            try
                Views = [Views, refl_Views];
            % If refl_Views empty because there are no valid views for this
            % number of reflections, then refl_Views will still be cleared
            % from last loop. This is okay so catch it and move on. If any
            % other error arises, rethrow it for the user to deal with.
            catch ME
                if ~strcmp(ME.identifier, 'MATLAB:refClearedVar')
                    rethrow(ME)
                end
            end
        end
        clear refl_Views
    end
    
%     assert(max_no_refl <= 2, "fn_make_geometry_views: for immersion, max number of geometry reflections must be <= 2.")
%     
%     view_el = 1;
%     for refl = 1:max_no_refl
%         
%         % Generate mode indices
%         char_modes = dec2bin([0:2^(refl+1)-1]);
%         view_modes = zeros(size(char_modes));
%         for col = 1:size(char_modes, 2)
%             view_modes(:, col) = str2num(char_modes(:, col)); %#ok<*ST2NM>
%         end
%         view_speeds1 = mat_speeds(view_modes+2);
%         view_speeds = couplant_spd * ones(size(view_speeds1, 1), size(view_speeds1, 2)+2);
%         view_speeds(:, 2:end-1) = view_speeds1;
% 
%         walls = 2*ones(refl, 1);
%         medium_ids = ones(refl+1, 1);
%         
%         % Generate wall indices
%         wall_idxs1 = dec2base([0:num_walls^refl-1], num_walls);
%         wall_idxs = zeros(size(wall_idxs1));
%         for col = 1:size(wall_idxs1, 2)
%             wall_idxs(:, col) = str2num(wall_idxs1(:, col)) + 1;
%         end
%         
%         % Get the valid ones.
%         wall_idxs1 = zeros(size(wall_idxs1));
%         this_row = 1;
%         for row = 1:size(wall_idxs1, 1)
%             is_adjacent = 0;
%             for col = 1:size(wall_idxs1, 2)-1
%                 if wall_idxs(row, col) == wall_idxs(row, col+1)
%                     is_adjacent = 1;
%                     break
%                 end
%             end
%             if ~is_adjacent
%                 wall_idxs1(this_row, :) = wall_idxs(row, :);
%                 this_row = this_row+1;
%             end
%         end
%         wall_idxs1(wall_idxs1(:, 1) == 0, :) = [];
%         fw_location_in_all_geom = find(where_F);
%         wall_idxs1(wall_idxs1 >= fw_location_in_all_geom) = wall_idxs1(wall_idxs1 >= fw_location_in_all_geom)+1;
%         wall_idxs = find(where_F) * ones(size(wall_idxs, 1), size(wall_idxs, 2)+2);
%         wall_idxs(:, 2:end-1) = wall_idxs1;
% 
%         view_geometries = all_geometries(wall_idxs);
%         
%         
% 
%         for wall = 1:size(view_geometries, 1)
%             for view = 1:size(view_modes, 1)
%                 % Get the name of this path.
%                 name = sprintf("%s", names(view_modes(view, 1)+1));
%                 rev_name = sprintf("%s", names(view_modes(view, end)+1));
%                 for mode = 2:size(view_modes, 2)
%                     name = sprintf("%s %s %s", name, view_geometries(wall, mode).name, names(view_modes(view, mode)+1));
%                     rev_name = sprintf("%s %s %s", rev_name, view_geometries(wall, end+1-mode).name, names(view_modes(view, end+1-mode)+1));
%                 end
%                 
%                 path_info = fn_path_info( ...
%                     name, ...
%                     rev_name, ...
%                     view_modes(view, :), ...
%                     view_geometries(wall, :).', ...
%                     view_speeds(view, :), ...
%                     mat_speeds, ...
%                     walls, ...
%                     medium_ids, ...
%                     densities, ...
%                     probe_freq, ...
%                     probe_pitch, ...
%                     probe_coords ...
%                 );
%                 
%                 if exist("Views", "var")==0
%                     Views = repmat(fn_create_geometry_view(path_info, all_geometries, probe_as_scatterer), tot_num_views, 1);
%                 else
%                     Views(view_el) = fn_create_geometry_view(path_info, all_geometries, probe_as_scatterer);
%                 end
%                 view_el = view_el + 1;
%             end
%         end
%     end
    
end

valid_view = zeros(size(Views));
for view = 1:size(Views, 1)
    valid_view(view) = any(Views(view).valid_path);
end

Views(~valid_view) = [];

% Get weights of valid views only.

view_walls = repmat("", size(Views));
view_modes = zeros(size(Views, 1), 2);

for view = 1:size(Views, 1)
%     Views(view).ray.weights = fn_compute_ray_weights(Views(view).ray, probe_freq);
    view_walls(view) = Views(view).ray.path_info.path_geometry.name;
    view_modes(view, :) = Views(view).ray.path_info.modes(:);
end

for view = 1:size(Views, 1)
    ray = Views(view).ray;
    if view_modes(view, 1) == view_modes(view, 2)
        inv_dir = ray.weights.directivity;
    else
        where_inv = logical(and(Views(view).ray.path_info.path_geometry.name == view_walls, ...
                            and(Views(view).ray.path_info.modes(1) == view_modes(:, 2), ...
                                Views(view).ray.path_info.modes(2) == view_modes(:, 1))));
        inv_dir = Views(where_inv).ray.weights.directivity;
    end
    [~, ~, num_freqs] = size(ray.weights.weights);
    Views(view).weights = zeros(probe_els^2, num_freqs);
    el = 1;
    for ii = 1 : probe_els
        for jj = 1 : probe_els
            Views(view).weights(el, :) = ( ...
                ray.weights.weights(ii, jj, :) * ray.weights.directivity(jj, ii, :) .* (1i * densities(medium_ids(end)+1) * view_speeds(view, end))...
            );

        el = el+1;
        end
    end
end
    
end
function ray = fn_compute_ray(scat_info, path_info, varargin)
% Function which computes the fermat path of a ray travelling from all
% probe positions to all scatterer positions. The Dijkstra method is used
% to do this, which is quicker than the bulk method when more than two legs
% are used.
%
% INPUTS:
% - scat_info : struct (1, 1)
%       Contains information on all scatterers which are being modelled.
%       Scatterer coordinates are in field `image_block`. Must be an output
%       from the fn_scat_info function.
% - path_info : struct(1, 1)
%       Contains information on the path which is being calculated. Must be
%       an output from the fn_path_info function.
% - geometry : OPTIONAL struct (no_walls, 1)
%       Contains all of the walls in the geometry which will be reflected
%       off. This is used for checking whether the rays are valid or not
%       due to collisions with walls.
% - freq_array : OPTIONAL array (no_freqs, 1)
%       Contains all of the frequencies which will be used to calculate ray
%       weights if the multi-frequency model is being used.
%
% OUTPUTS:
% - ray : struct (1, 1)
%       Contains all information calculated from the ray, including the
%       path taken and the ray weights.

% Unpack scat_info and path_info structures.
scatterers = scat_info.image_block;
path_geometry = path_info.path_geometry;
speeds = path_info.speeds;
probe_coords = path_info.probe_coords;
[probe_els, ~] = size(probe_coords);
[num_scatterers, ~] = size(scatterers);

% Work out the number of walls, depending on whether we are in the contact
% or immersion case.
if ~isstruct(path_geometry) % We are in the contact case.
    no_walls = 0;
else % We are in the immersion case.
    [no_walls, ~] = size(path_geometry);
    [wall_pixels, ~] = size(path_geometry(1).coords);
end

% Initialise the structure.
ray.min_times = zeros(probe_els, num_scatterers);
ray.wall_idxs = zeros(probe_els, num_scatterers, no_walls);
ray_coords = zeros(probe_els, num_scatterers, no_walls+2, 3);
ray.path_info = path_info;
ray.scat_info = scat_info;
ray.valid_paths = zeros(probe_els, 1);

% Calculate the fermat path.
for scat = 1 : num_scatterers
    for tx = 1 : probe_els
        % In the contact case, we trace directly from probe to scatterer.
        % Check whether we are dealing with contact or immersion.
        if ~isstruct(path_geometry) % If we are in the direct contact case.
            
            % There is only one leg. Fermat path == straight line from
            % probe to scatterer.
            min_times = sqrt( ...
                (scatterers(scat, 1) - probe_coords(tx, 1)) ^ 2 + ...
                (scatterers(scat, 3) - probe_coords(tx, 3)) ^ 2) / ...
            speeds(1);
            
            % Get the ray coordinates for use with fn_valid_paths.
            ray_coords(tx, scat, 1, :) = probe_coords(tx, :);
            ray_coords(tx, scat, 2, :) = scatterers(scat, :);
            
            % Collect the values.
            ray.min_times(tx, scat) = min_times;
        
            
        % We must be in the skip contact case, or the immersion case. In
        % other words, we will have more than one leg, and thus we need to
        % use the Dijkstra method to calculate the fermat path.
        else 
            
            % Initialise the arrays for minimum times and associated wall
            % indices (i.e. the index of the discrete point on the wall
            % which the ray is transmitted through/reflected from.
            min_times = zeros(wall_pixels, 1);
            wall_idxs = linspace(1, wall_pixels, wall_pixels);

            % First and last legs of the ray are simple to calculate. Treat
            % them differently. Start with first leg.
            
            % Calculate the time taken to travel from the probe to all
            % points on the first wall.
            min_times(:, 1) = ( ...
                sqrt((path_geometry(1).coords(:, 1) - probe_coords(tx, 1)) .^ 2 + ...
                     (path_geometry(1).coords(:, 3) - probe_coords(tx, 3)) .^ 2) ./ speeds(1) ...
            );

            % If more than one boundary is present, then we will have
            % middle legs that we must consider - this is where the
            % efficiency improvement of Dijkstra comes in.
            if no_walls > 1
                % Compute the times from all points on the previous wall to
                % all points on the next wall, i.e. compute a total of 
                % (nth_wall_points)*((n+1)th_wall_points) times. Add these
                % to the corresponding previously calculated times (i.e.
                % the time from the probe to the kth point on the nth wall,
                % added to the time calculated from the kth point on the
                % nth wall to all points on the (n+1)th wall.
                for wall = 2 : no_walls
                    min_times_matrix = min_times + ...
                        sqrt((path_geometry(wall).coords(:, 1).' - path_geometry(wall-1).coords(:, 1)).^2 + ...
                             (path_geometry(wall).coords(:, 3).' - path_geometry(wall-1).coords(:, 3)).^2) ./ speeds(wall);
                    
                % Find the minimum time from the probe to each point on the
                % (n+1)th wall, and return these times. This min_times
                % object is now equivalent to the min_times object before
                % this leg, except it has now been propagated one wall
                % further.
                    [find_i, find_j] = ind2sub(size(min_times_matrix), find(min_times_matrix == min(min_times_matrix, [], 2)));
                    for k = 1 : wall_pixels
                        min_times(k, 1) = min_times_matrix(find_i(k), find_j(k));
                    end
                    idxs = wall_idxs(:, find_i);
                    wall_idxs = zeros(wall, wall_pixels);
                    wall_idxs(1:wall-1, :) = idxs;
                    wall_idxs(wall, :) = linspace(1, wall_pixels, wall_pixels);
                end
            end

            % Finally, compute the last leg.
            min_times = min_times + ...
                sqrt((scatterers(scat, 1) - path_geometry(end).coords(:, 1)) .^ 2 + ...
                     (scatterers(scat, 3) - path_geometry(end).coords(:, 3)) .^ 2) ./ speeds(no_walls + 1);
            
            % Now that we have reached the scatterer, find the path which
            % takes the least time. Record this time, and the wall indices
            % which describe the path taken through the geometry.
            find_i = find(min_times == min(min_times, [], 'all'));
            if size(find_i, 1) > 1
                find_i = find_i(1);
            end
            ray.min_times(tx, scat) = min_times(find_i);
            ray.wall_idxs(tx, scat, 1:no_walls) = wall_idxs(:, find_i);
            
            ray_coords(tx, scat, 1, :) = probe_coords(tx, :);
            for wall = 1:no_walls
                ray_coords(tx, scat, wall+1, :) = path_geometry(wall).coords(ray.wall_idxs(tx, scat, wall), :);
            end
            ray_coords(tx, scat, end, :) = scatterers(scat, :);
                
        end
    end
end



if nargin > 2
    for arg = 1 : nargin-2
        argument = varargin{arg};
        if isa(argument, 'double')
            freq_array = varargin{arg};
            % If frequency is provided, compute the ray weights.
            ray.weights = fn_compute_ray_weights(ray, freq_array);
            ray.freq_array = freq_array;
        elseif isa(argument, 'struct')
            geometry = varargin{arg};
            % Determine whether the ray paths are valid.
            ray.valid_paths = fn_valid_paths(path_info, ray_coords, geometry);
        end
    end
end

end
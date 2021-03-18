function valid_paths = fn_valid_paths(path_info, ray_coords, all_geometry)
% Determines whether the path traced by each ray is valid by detecting
% whether the legs of each ray intersect with any of the walls in the
% geometry. Note that here, "intersect" means that there is at least one
% point which is coincident between the line representing the leg of a ray
% and a wall, not including the points where reflections occur.
%
% Note that this function also checks whether each path passes through an
% extreme point on a wall (i.e. a corner), and treats the path as invalid
% if it does. These paths will always be unphysical, or if they technically
% are valid then they will be governed by different physics than are
% currently accounted for.
%
% INPUTS:
%   - path_info : struct (1, 1)
%       Contains information on the path which is being calculated. Must be
%       an output from the fn_path_info function.
%   - ray_coords : array (probe_els, num_scatterers, num_walls+2, 3)
%       The coordinates of the ray at the start and end of each leg. As the
%       end of one leg has the same coordinates as the start of the next,
%       there will be num_legs+1 coordinates describing the ray, or
%       num_walls+2.
%   - all_geometry : struct (total_num_walls, 1)
%       All walls which have been specified to define the inspection block.
%       This must at least contain the walls from which the ray reflects,
%       but may contain extra ones which the ray does not interact with in
%       this path.
%
% OUTPUTS:
%   - valid_paths : array (probe_els, num_scatterers)
%       Array of logical values which indicate whether the ray path
%       (described by wall_idxs) is valid or not, i.e. other than the walls
%       which are reflected from / transmitted through, does the ray
%       intersect any other walls?

% Get useful dimensions.
[probe_els, num_scatterers, num_points, ~] = size(ray_coords);
num_legs = num_points - 1;

% Get the geometry which the wall reflects from. We'll be checking for
% intersections with walls, but the ray will already intersect walls from
% which it transmits/reflects. To exclude these walls from intersection
% checks, get the ones we need to exclude here.
path_geometry = path_info.path_geometry;

% Initialise output array. Assume valid to start, and then change if
% intersection is found.
valid_paths = ones(probe_els, num_scatterers);

% Each ray traces from a probe element to a scatterer.
for el = 1:probe_els
    for scat = 1:num_scatterers
        % For each leg in the current ray.
        for leg = 1:num_legs
            % Get the start and end point of each leg.
            leg_start = squeeze(ray_coords(el, scat, leg, :));
            leg_end = squeeze(ray_coords(el, scat, leg+1, :));
            
            % Get the walls we need to check for intersections with. We may
            % need to exclude a certain number of walls depending on the
            % path and which leg of the path we are currently on.
            %
            % If we are on a direct contact path, we do not reflect from
            % any geometry, thus we test for intersections with all
            % geometry.
            if num_legs == 1
                geometry_for_testing = all_geometry;
                
            % If we are not (i.e. more than 1 leg in total) then check if
            % we are on the first or the last leg. If we are, there is only
            % one wall to exclude.
            elseif leg == 1 % If we are on the first leg.
                % First check whether we pass through a corner. On the
                % first leg, leg_start is on the probe and leg_end is on
                % the wall.
                if or(leg_end.' == path_geometry(leg).point1, leg_end.' == path_geometry(leg).point2)
                    valid_paths(el, scat) = 0;
                    break
                end
                
                % If the path is unphysical, we break out of the loop over
                % the path and never reach this point. If it is physical,
                % then carry on and check for intersections.
                geometry_for_testing = all_geometry;
                exclude_name = path_geometry(leg).name;
                for wall = 1:size(all_geometry, 1)
                    if geometry_for_testing(wall).name == exclude_name
                        geometry_for_testing(wall) = [];
                        break
                    end
                end
                
            elseif leg == num_legs % If we are on the last leg.
                % Check if the path goes through a wall corner.
                if or(leg_start.' == path_geometry(leg-1).point1, leg_end.' == path_geometry(leg-1).point2)
                    valid_paths(el, scat) = 0;
                    break
                end
                
                % Exclude the wall we interact with.
                geometry_for_testing = all_geometry;
                exclude_name = path_geometry(leg-1).name;
                for wall = 1:size(all_geometry, 1)
                    if geometry_for_testing(wall).name == exclude_name
                        geometry_for_testing(wall) = [];
                        break
                    end
                end
                
            % We must be on a path with >= 3 legs, and the leg we are
            % currently looking at must not be an end leg.
            else
                % Check if the path goes through a wall corner.
                if or(leg_end.' == path_geometry(leg).point1, leg_end.' == path_geometry(leg).point2)
                    valid_paths(el, scat) = 0;
                    break
                elseif or(leg_start.' == path_geometry(leg-1).point1, leg_start.' == path_geometry(leg-1).point2)
                    valid_paths(el, scat) = 0;
                    break
                end
                
                % Exclude the wall we interact with.
                geometry_for_testing = all_geometry;
                exclude_name = path_geometry(leg-1).name;
                for wall = 1:size(all_geometry, 1)
                    if geometry_for_testing(wall).name == exclude_name
                        geometry_for_testing(wall) = [];
                        break
                    end
                end
                exclude_name = path_geometry(leg).name;
                for wall = 1:size(all_geometry, 1)
                    if geometry_for_testing(wall).name == exclude_name
                        geometry_for_testing(wall) = [];
                        break
                    end
                end
            end
            
            % Now that we have the geometry to test for intersections with,
            % we can work out whether this path is valid.
            %
            % Check that there is some geometry remaining to test with. If 
            % there is no geometry, then the path will automatically be
            % valid, so move on.
            if size(geometry_for_testing, 2) ~= 0
                % Check that all walls do not intersect.
                for wall = 1:size(geometry_for_testing, 1)
                    % If there is an intersection, set the validity to 0
                    % and break out of the loop over this ray's legs, as we
                    % do not need to check any more.
                    wall_start = geometry_for_testing(wall).coords(1, :);
                    wall_end = geometry_for_testing(wall).coords(end, :);
                    is_intersection = fn_is_intersection(wall_start, wall_end, leg_start, leg_end);
                    if is_intersection
                        valid_paths(el, scat) = 0;
                        break
                    end
                end
                % Break out of the loop over this ray's legs.
                if is_intersection
                    break
                end
            
            % If there are no walls to test, then there will be no
            % intersections as all walls are straight lines. Move on to the
            % next path.
            else
                break
            end
            
        end
    end
end



end
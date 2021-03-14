function valid_paths = fn_valid_paths(path_info, ray_coords, all_geometry)
% Determines whether the path traced by each ray is valid by detecting
% whether the legs of each ray intersect with any of the walls in the
% geometry. Note that here, "intersect" means that there is at least one
% point which is coincident between the line representing the leg of a ray
% and a wall, not including the points where reflections occur.
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

% Get path geometry so that we know which walls we need to consider
% intersections with (i.e. we want the subset of all_geometry which
% excludes the walls that the specific leg we are dealing with is adjoined
% by.
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
            
            % Get the walls we need to exclude from all_geometry so that we
            % can test for intersection. The way we exclude geometry will
            % depend on the total number of legs, as well as whether we are
            % in the first or last leg.
            % 
            % If we are on a direct contact path, we do not reflect from
            % any geometry. Test for intersection against all_geometry.
            if num_legs == 1
                geometry_for_testing = all_geometry;
                
            % If we are not (i.e. more than 1 leg in total) then check if
            % we are on the first or the last leg. If we are, there is only
            % one wall to exclude.
            elseif leg == 1 % If we are on the first leg.
                geometry_for_testing = all_geometry;
                exclude_name = path_geometry(leg).name;
                for wall = 1:size(all_geometry, 1)
                    if geometry_for_testing(wall).name == exclude_name
                        geometry_for_testing(wall) = [];
                        break
                    end
                end
                
            elseif leg == num_legs % If we are on the last leg.
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
            % we can work out whether this is valid or not.
            %
            % Check that there is some geometry remaining to test with. Ray
            % must be valid if not, so move on.
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
            
            % If there are no walls to test, then there can be no
            % intersections. This comes from the fact that walls must be
            % straight lines, and thus we must either have no walls at all,
            % or we are reflecting off all the available walls. If we are
            % reflecting, then the ray cannot intersect with this wall.
            else
                break
            end
            
        end
    end
end



end
function geometry = fn_make_geometry(is_contiguous, profile_list, N, varargin)
% Assembles all of the walls in the geometry, from which rays will be
% traced. Note that walls may be considered contiguous (e.g. C-shaped
% geometry with front, back and side wall for use in immersion), or not
% (e.g. shaped like an =, with front and back wall only and open ends).
%
% INPUTS:
% - is_contiguous : logical
%       Logical test for whether points provided are contiguous. When
%       contiguous, geometry does not have to be closed (e.g. ☐-shape and
%       C-shape both allowed).
% - profile_list : array of strings (M, 1), or 0.
%       List of wall profiles, must be equal to the number of walls given
%       (i.e. if is_contiguous==true then length must be equal to number of
%       points - 1; if is_contiguous==false then length must be equal to
%       the number of points / 2). Valid profiles are "line", "arcL" for an
%       arc of a circle which curves left, and "arcR" for an arc which
%       curves right. If wall_profile_list==0, then all walls are assumed
%       to be lines.
% - N : integer
%       Number of points each wall will be discretised into. This is a
%       function of final TFM pixel size, and should be fine enough that
%       rays traced from this wall are approximately the same as that
%       predicted by Snell's law. Needs validation, but keep wall pixel
%       size ~= 0.1 * TFM pixel size.
% - points : double array (M, 3)
%       Either an array of shape specified listing all of the points in the
%       geometry, or each point is specified as an additional argument e.g.
%       fn_make_geometry(..., [1,0,0], [1,0,-1], [-1,0,-1], [-1,0,0])
%
% OUTPUTS:
% - geometry : struct (N, 1)
%       Set of all walls in the geometry. Each element is output from the
%       fn_make_wall function.

% Assembles all of the walls in the geometry, from which rays will be
% reflected. Note that this function can assume that the walls are
% contiguous (e.g. C-shaped geometry with front, back and sidewall) or not
% (e.g. =-shaped geometry with front and backwall only).

VALID_PROFILES = ["line", "arcL", "arcR"];

if and(is_contiguous ~= 0, is_contiguous ~= 1)
    error("fn_make_geometry: is_contiguous not logical")
end
% Are there enough arguments?
if nargin < 4
    error("fn_make_geometry: insufficient number of points")
end

% If points provided as single array
if nargin == 4
    % Is there at least one wall?
    if size(varargin{1}, 1) < 2
        error("fn_make_geometry: insufficient number of points")
    end
    % If not contiguous, do points come in pairs?
    if and(~is_contiguous, ~mod(size(varargin{1}, 1), 2))
        error("fn_make_geometry: insufficient number of points")
    end
    
    corners = varargin{1};
    
% If points provided as separate arguments
else
    % Is there at least one wall?
    if nargin-3 < 2
        error("fn_make_geometry: insufficient number of points")
    end
    % If not contiguous, do points come in pairs?
    if and(~is_contiguous, ~mod(nargin-3, 2))
        error("fn_make_geometry: insufficient number of points")
    end
    
    corners = zeros(nargin-3, 3);
    for point = 1:nargin-3
        coords = varargin{point};
        if size(coords, 2) ~= 3
            error("fn_make_geometry: points have wrong number of dimensions")
        end
        corners(point, :) = coords;
    end
end

% Generate profiles if not given.
if isa(profile_list, 'double')
    profile_list = repmat("line", size(corners, 1), 1);
% If profiles provided, are there enough?
else
    if is_contiguous
        if size(profile_list, 2) ~= size(corners, 1)-1
            error("fn_make_geometry: not enough profiles specified")
        end
    else
        if size(profile_list, 1) ~= size(corners, 1)/2
            error("fn_make_geometry: not enough profiles specified")
        end
    end
    % Check that all profiles are valid.
    for profile = profile_list
        if ~any(strcmpi(profile, VALID_PROFILES))
            error("fn_make_geometry: invalid profile %s", profile)
        end
    end
end

side_count   = 0;
back_count   = 0;
other_count  = 0;
curved_count = 0;

% If the walls are not contiguous, treat each pair of points as the ends of
% the wall.
if ~is_contiguous
    no_walls = size(corners, 1)/2;
    
    % Pre-generate geometry to overwrite later.
    geometry = repmat(fn_make_wall('name', [0,0,0], [0,0,0], 1, 1), no_walls, 1);
    
    for wall = 1:no_walls
        start_point = corners(2*wall - 1, :);
        end_point = corners(2*wall, :);
        
        wall_id = 0;
        
        if strcmpi(profile_list(wall), "line")
            if and(start_point(3) == -1e-5, end_point(3) == -1e-5)
                name = sprintf("%s", "F");
                wall_id = 1;
            elseif (start_point(3) - end_point(3)) == 0
                back_count = back_count + 1;
                name = sprintf("%s%d", "B", back_count);
                wall_id = 2;
            elseif (start_point(1) - end_point(1)) == 0
                side_count = side_count + 1;
                name = sprintf("%s%d", "S", side_count);
                wall_id = 3;
            else
                other_count = other_count + 1;
                name = sprintf("%s%d", "O", other_count);
            end
        
            geometry(wall) = fn_make_line(name, start_point, end_point, N, wall_id);
        else
            curved_count = curved_count + 1;
            name = sprintf("%s%d", "C", curved_count);
            profile = char(profile_list(wall));
            if strcmpi(profile(4), "L")
                go_left = true;
            elseif strcmpi(profile(4), "R")
                go_left = false;
            end
            geometry(wall) = fn_make_arc(name, start_point, end_point, 15e-3, N, go_left, wall_id); 
        end
    end
    
% If the walls are contiguous.
else
    no_walls = size(corners, 1) - 1;

    geometry = repmat(fn_make_line('name', [0,0,0], [0,0,1], 1, 1), no_walls, 1);
    
    for wall = 1:no_walls
        start_point = corners(wall, :);
        end_point = corners(wall + 1, :);
        
        wall_id = 0;
        
        if strcmpi(profile_list(wall), "line")
            if and(start_point(3) == -1e-5, end_point(3) == -1e-5)
                name = sprintf("%s", "F");
                wall_id = 1;
            elseif (start_point(3) - end_point(3)) == 0
                back_count = back_count + 1;
                name = sprintf("%s%d", "B", back_count);
                wall_id = 2;
            elseif (start_point(1) - end_point(1)) == 0
                side_count = side_count + 1;
                name = sprintf("%s%d", "S", side_count);
                wall_id = 3;
            else
                other_count = other_count + 1;
                name = sprintf("%s%d", "O", other_count);
            end
            
            geometry(wall) = fn_make_line(name, start_point, end_point, N, wall_id);
        else
            curved_count = curved_count + 1;
            name = sprintf("%s%d", "C", curved_count);
            profile = char(profile_list(wall));
            if strcmpi(profile(4), "L")
                go_left = true;
            elseif strcmpi(profile(4), "R")
                go_left = false;
            end
            geometry(wall) = fn_make_arc(name, start_point, end_point, 15e-3, N, go_left, wall_id);      
        end
    end
end

wall_names = repmat("", size(geometry, 1), 1);
for wall = 1:size(geometry, 1)
    wall_names(wall) = geometry(wall).name;
end
where_F = logical(count(wall_names, "F"));
assert(sum(where_F)<=1, "fn_make_geometry: too many frontwalls");

end
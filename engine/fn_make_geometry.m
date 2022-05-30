function geometry = fn_make_geometry(is_contiguous, N, varargin)
% Assembles all of the walls in the geometry, from which rays will be
% reflected. Note that this function can assume that the walls are
% contiguous (e.g. C-shaped geometry with front, back and sidewall) or not
% (e.g. =-shaped geometry with front and backwall only).
%
% INPUTS:
% - is_contiguous : logical
%       Logical test of whether the provided points will come together to
%       form a contiguous geometry or not. Note that when contiguous, the
%       geometry does not have to be closed (e.g. ☐-shaped), and can be
%       open (C-shaped).
% - N : integer
%       Number of points to discretise walls into. Note that this should be
%       a function of the expected final TFM pixel size: if wall pixels are
%       >= image pixels, then it is very likely that there will be extreme
%      inaccuracy in rays which reflect close to the wall. Try to keep wall
%      pixels at about 0.1 * image pixels (not checked).
% - point : array (n, 3)
%       3D coordinates of the ends of walls. When is_contiguous == 0, an
%       even number of points must be provided. Each pair of points are
%       then treated as the ends of the walls. When is_contiguous == 1, N+1
%       points must be provided, where N is the number of walls to be
%       modelled. Each point is then treated as a corner of the geometry.
%
% OUTPUTS:
% - geometry : struct (N, 1)
%       Set of all walls in the geometry. Each element is output from the
%       fn_make_wall function.

assert(or(is_contiguous==0, is_contiguous==1), ...
    "fn_make_geometry: is_contiguous is not logical.")
assert(nargin >= 3, ...
    "fn_make_geometry: Insufficient number of points.")
corners = zeros(0, 3);
for point = 1 : nargin-2
    coords = varargin{point};
    assert(size(coords, 2) == 3, ...
        "fn_make_geometry: Points have invalid number of dimensions.")
    corners = [corners; coords];
end

side_count = 0;
back_count = 0;
other_count = 0;

% If the walls are not contiguous, treat each pair of points as the ends of
% the wall.
if ~is_contiguous
    assert(int8(size(corners, 1)/2) == size(corners, 1)/2, ...
        "fn_make_geometry: Incorrect number of points.")
    no_walls = size(corners, 1)/2;

    geometry = repmat(fn_make_wall('name', [0,0,0], [0,0,0], 1, 1), no_walls, 1);
    
    for wall = 1:no_walls
        start_point = corners(2*wall - 1, :);
        end_point = corners(2*wall, :);
%         start_point = varargin{2*wall - 1};
%         end_point = varargin{2*wall};
        
        wall_id = 0;
        
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
        
        geometry(wall) = fn_make_wall(name, start_point, end_point, N, wall_id);
    end
    
% If the walls are contiguous.
else
    no_walls = size(corners, 1) - 1;

    geometry = repmat(fn_make_wall('name', [0,0,0], [0,0,0], 1, 1), no_walls, 1);
    
    for wall = 1:no_walls
        start_point = corners(wall, :);
        end_point = corners(wall + 1, :);
%         start_point = varargin{wall};
%         end_point = varargin{wall + 1};
        
        wall_id = 0;
        
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
        
        geometry(wall) = fn_make_wall(name, start_point, end_point, N, wall_id);
    end
end

wall_names = repmat("", size(geometry, 1), 1);
for wall = 1:size(geometry, 1)
    wall_names(wall) = geometry(wall).name;
end
where_F = logical(count(wall_names, "F"));
assert(sum(where_F)<=1, "fn_make_geometry: too many frontwalls");

end
function geometry = fn_make_line(name, point1, point2, N, wall_id)
% Creates a single wall object for simulating TFMs and sensitivity maps.
% Currently only supports planar walls in 2D, and walls are assumed to be
% in the x-z plane (i.e. y == 0).
%
% INPUTS:
% - name : string
%       Name of the wall.
% - point1 : array (1, 3)
%       3D coordinates of a point which defines one end of the wall.
% - point2 : array (1, 3)
%       3D coordinates of a point which defines one end of the wall.
% - N : integer
%       The number of points which the wall will be discretised over.
%
% OUTPUTS:
% - geometry : struct
%      Structure of shape (1, 1) which contains all defining information of
%      a single wall in the geometry.

assert(point1(2) == 0, "Wall is not located in the x-z plane.")
assert(point2(2) == 0, "Wall is not located in the x-z plane.")

% Assign point1 and point2 correctly so that the bases are always in the
% same direction for the same wall.
points = [point1; point2];
[~, I] = sort(points(:, 1), 'descend');
points = points(I, :);
[~, I] = sort(points(:, 3), 'ascend');
points = points(I, :);

point1 = points(1, :);
point2 = points(2, :);

geometry.name = name;
geometry.point1 = point1;
geometry.point2 = point2;

geometry.coords = zeros(N+1, 3);
geometry.coords(:, 1) = linspace(point1(1), point2(1), N+1);
geometry.coords(:, 2) = linspace(point1(2), point2(2), N+1);
geometry.coords(:, 3) = linspace(point1(3), point2(3), N+1);
geometry.coords = geometry.coords(1:end-1, :);

% The basis is the vector which define the rotation of the geometry.
% basis(:, 1) is parallel to the wall; basis(:, 2) should always be
% parallel to the y-axis; basis(:, 3) is normal to the wall. This normal
% vector is used to work out the angle each ray makes with the wall.
%
% In arim, when multiple (parallel) walls are considered, the 3rd basis
% vector should point in the same direction for all walls. Now that we're
% considering multiple angles of wall geometries, this needs examination.
% For now, have 3rd basis vector (i.e. the one normal to the surface)
% pointing in the +ve z, +ve x direction (in the standard basis), with z
% taking priority.
geometry.basis = zeros(1, 3, 3);
geometry.basis(:, :, 1) = (point1 - point2) / norm(point2 - point1);
geometry.basis(:, :, 2) = [0, 1, 0];
geometry.basis(:, :, 3) = [-geometry.basis(:, 3, 1), geometry.basis(:, 2, 1), geometry.basis(:, 1, 1)];

geometry.wall_id = wall_id;

end
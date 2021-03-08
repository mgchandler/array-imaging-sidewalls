function geometry = fn_make_wall(name, point1, point2, N, wall_id)
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
geometry.basis = zeros(3, 3);
geometry.basis(:, 1) = (point1 - point2) / norm(point1 - point2);
geometry.basis(:, 2) = [0, 1, 0];
geometry.basis(:, 3) = [-geometry.basis(3, 1), geometry.basis(2, 1), geometry.basis(1, 1)];

geometry.wall_id = wall_id;

end
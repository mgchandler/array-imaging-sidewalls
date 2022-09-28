function geometry = fn_make_arc(name, point1, point2, rad, N, go_left, wall_id)
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
% - rad : double
%       Radius of the arc.
% - N : integer
%       The number of points which the wall will be discretised over.
% - go_left : logical
%       Along the line from point1 -> point2, does the arc curve to the
%       left?
%
% OUTPUTS:
% - geometry : struct
%      Structure of shape (1, 1) which contains all defining information of
%      a single wall in the geometry.

assert(point1(2) == 0, "Wall is not located in the x-z plane.")
assert(point2(2) == 0, "Wall is not located in the x-z plane.")

% % Assign point1 and point2 correctly so that the bases are always in the
% % same direction for the same wall.
% points = [point1; point2];
% [~, I] = sort(points(:, 1), 'descend');
% points = points(I, :);
% [~, I] = sort(points(:, 3), 'ascend');
% points = points(I, :);
% 
% point1 = points(1, :);
% point2 = points(2, :);

geometry.name = name;
geometry.point1 = point1;
geometry.point2 = point2;

% point1 and point2 sit on the arc by definition. Find centre point using
% radius.
half_chord = (point2 - point1) / 2;
normal = (-1)^go_left * sqrt(rad^2 - norm(half_chord)^2) * [-half_chord(3)/norm(half_chord), 0, half_chord(1)/norm(half_chord)];
% Normal vector defined wrt half_chord vector; which is defined wrt point1.
centre = point1 + half_chord + normal;

% Angles swept through from point1 -> point2.
phi1 = atan2(point1(3) - centre(3), point1(1) - centre(1));
phi2 = atan2(point2(3) - centre(3), point2(1) - centre(1));
phi  = linspace(phi1, phi2, N+1);

geometry.coords = zeros(N+1, 3);
geometry.coords(:, 1) = rad * real(exp(1i*phi)) + centre(1);
geometry.coords(:, 2) = zeros(N+1, 1);
geometry.coords(:, 3) = rad * imag(exp(1i*phi)) + centre(3);
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
geometry.basis = zeros(N+1, 3, 3);
geometry.basis(:, :, 1) = [cos(phi).', zeros(N+1,1), sin(phi).'];
geometry.basis(:, :, 2) = repmat([0,1,0], N+1, 1);
geometry.basis(:, :, 3) = [-geometry.basis(:, 3, 1), zeros(N+1, 1), geometry.basis(:, 1, 1)];

geometry.wall_id = wall_id;

end
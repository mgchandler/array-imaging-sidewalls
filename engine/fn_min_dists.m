function dists = fn_min_dists(geometries)
% Calculates the distances of each leg of the ray, returning x, y, z and
% Euclidean distance of each leg. Only returns one set of distances,
% assuming that loops over probe and scatterer are in the parent function.
% N.B. x, y and z can be -ve if start -> end runs in the -ve direction, but
% Euclid. dist. will always be +ve.
%
% INPUTS:
% - geometries : array of size (no_walls+2, dims)
%       List of coordinates on the walls through which the ray passes from
%       probe to scatterer. Note that the first point must be the probe, 
%       the middle points must be the geometries (in order) and the final 
%       point must be the scatterer.
%
% OUTPUTS:
% - dists : array (no_legs, 4)
%       The distances of each leg in the provided ray. In the 2nd
%       dimension, the distances in the x, y and z axes are provided, as
%       well as the Euclidean distance.

[no_geoms, ~] = size(geometries);
no_legs = no_geoms-1;

dists = zeros(no_legs, 4);

for leg = 1:no_legs
    dists(leg, 1) = geometries(leg+1, 1) - geometries(leg, 1);
    dists(leg, 2) = geometries(leg+1, 2) - geometries(leg, 2);
    dists(leg, 3) = geometries(leg+1, 3) - geometries(leg, 3);
    dists(leg, 4) = sqrt(dists(leg, 1)^2 + dists(leg, 2)^2 + dists(leg, 3)^2);
end

end
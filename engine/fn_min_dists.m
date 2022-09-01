function dists = fn_min_dists(geometries)
% Calculates the distances of each leg of the ray, returning x, y, z and
% Euclidean distance of each leg. Only returns one set of distances,
% assuming that loops over probe and scatterer are in the parent function.
% N.B. x, y and z can be -ve if start -> end runs in the -ve direction, but
% Euclid. dist. will always be +ve.
%
% INPUTS:
% - geometries : array of size (no_walls+2, dims, no_scats)
%       List of coordinates on the walls through which the ray passes from
%       probe to scatterer. Note that the first point must be the probe, 
%       the middle points must be the geometries (in order) and the final 
%       point must be the scatterer.
%
% OUTPUTS:
% - dists : array (no_legs, dims+1, no_scats)
%       The distances of each leg in the provided ray. In the 2nd
%       dimension, the distances in the x, y and z axes are provided, as
%       well as the Euclidean distance.

[no_geoms, dims, no_scats] = size(geometries);
no_legs = no_geoms-1;

dists = zeros(no_legs, dims+1, no_scats);

for leg = 1:no_legs
    for dim = 1:dims
        dists(leg, dim, :) = geometries(leg+1, dim, :) - geometries(leg, dim, :);
    end
    dists(leg, end, :) = sqrt(sum(dists(leg, :, :).^2, 2));
end

end
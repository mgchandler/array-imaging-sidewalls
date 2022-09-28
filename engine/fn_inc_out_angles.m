function alpha_beta = fn_inc_out_angles(dists, path_geometry, wall_idxs)
% Computes the incident angles (alpha) and outgoing angles (beta) on each
% wall. Assumes that the loop for probes and scatterers is in the parent
% function.

% INPUTS:
% - dists : array of size (no_legs, dims+1, no_scats)
%       Distances between each wall.
% - path_geometry : struct (no_walls, 1)
%       Set of geometry which the ray will interact with. Does not include
%       the probe or the scatterer.
% - wall_idxs : array (1, no_walls)
%       Indices of the wall point which the ray from element tx passes
%       through. Equivalent to wall_idxs field in output of fn_compute_ray.
%
% OUTPUTS:
% - alpha_beta : array (no_legs, 4, no_scats)
%       The incident and outgoing angles made for each leg. In the 2nd
%       axis, the order of angles provided is:
%           [conv_inc_ang, conv_out_ang, sign_inc_ang, sign_out_ang]
%       where 'conv' is conventional, and 'sign' is signed, both names
%       co-opted from `arim`. In the 1st axis, the angles provided are the
%       incident angles on each wall and the scatterer, and the outgoing
%       angles from the probe and walls.

[no_legs, dims_plus_one, no_scats] = size(dists);
alpha_beta = zeros(no_legs, 4, no_scats);
for wall = 1:no_legs-1
    % Flat walls only have one basis.
    if size(path_geometry(wall).basis, 1) == 1
        wall_idxs = 1;
    % If the wall is curved, make sure that idxs are valid.
    else
        if wall_idxs(wall) > size(path_geometry(wall).basis, 1)
            error("fn_inc_out_angles: idx ray passes through too large for wall.")
        end
    end
end

% Compute incident angles on walls and scatterer.
for leg = 1:no_legs
    % Outgoing angle on 1st leg and incident angle on last leg have no
    % corresponding path geometry to reference. Treat these differently.
    
    %% Incident angles.
    
    % Get the vector representing this leg in the wall's coordinate system.
    %
    % In arim, angles are worked out from the coordinates of the source
    % point of the ray (or the end of the ray in the outgoing case).
    % Here, dists is the direction of the ray. To get the source
    % points, do source = current - ray_direction. To get the angle, we
    % can treat `current` as the origin, so simply take the -ve of the
    % ray direction.
    if leg ~= no_legs
        dist = reshape(-dists(leg, 1:dims_plus_one-1, :), dims_plus_one-1, 1, no_scats);
        cart_inc = zeros(dims_plus_one-1, no_scats);
        for scat = 1:no_scats
            cart_inc(:, scat) = linsolve(squeeze(path_geometry(leg).basis(wall_idxs(leg), :, :)), dist(:, :, scat));
        end
    else
        cart_inc = reshape(-dists(leg, 1:dims_plus_one-1, :), dims_plus_one-1, no_scats);
    end
    
    rad_inc = reshape(dists(leg, dims_plus_one, :), 1, no_scats);
    % Compute angle.
    alpha = acos(cart_inc(3, :) ./ rad_inc);
    
    % Work out how to shift signed angle.
    azimuth = atan2(cart_inc(2, :), cart_inc(1, :));
    sign_test = and(-pi/2 < azimuth, azimuth <= alpha);
    sign_alpha = alpha .* (2*sign_test-1);
    
    % Work out how to shift conventional ang
    conv_test = cart_inc(3, :) < 0;
    conv_alpha = conv_test*pi + (-1).^conv_test .* alpha;
    
    alpha_beta(leg, 1, :) = conv_alpha;
    alpha_beta(leg, 3, :) = sign_alpha;
    
    %% Outgoing angles
    
    % Get the vector representing this leg in the wall's coordinate system.
    if leg ~= 1
        dist = reshape(-dists(leg, 1:dims_plus_one-1, :), dims_plus_one-1, 1, no_scats);
        cart_inc = zeros(dims_plus_one-1, no_scats);
        for scat = 1:no_scats
            cart_inc(:, scat) = linsolve(squeeze(path_geometry(leg-1).basis(wall_idxs(leg-1), :, :)), dist(:, :, scat));
        end
    else
        cart_out = reshape(-dists(leg, 1:dims_plus_one-1, :), dims_plus_one-1, no_scats);
    end
    
    rad_out = reshape(dists(leg, dims_plus_one, :), 1, no_scats);
    % Compute angle.
    beta = acos(cart_out(3, :) ./ rad_out);
    % Work out how to shift signed angle.
    azimuth = atan2(cart_out(2, :), cart_out(1, :));
    sign_test = and(-pi/2 < azimuth, azimuth <= beta);
    sign_beta = beta .* (2*sign_test-1);
    % Work out how to shift conventional ang
    conv_test = cart_inc(3, :) < 0;
    conv_beta = conv_test*pi + (-1).^conv_test .* beta;
    
    alpha_beta(leg, 2, :) = conv_beta;
    alpha_beta(leg, 4, :) = sign_beta;
    
end

end
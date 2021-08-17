function alpha_beta = fn_inc_out_angles(dists, path_geometry)
% Computes the incident angles (alpha) and outgoing angles (beta) on each
% wall. Assumes that the loop for probes and scatterers is in the parent
% function.

% INPUTS:
% - dists : array of size (no_legs, 4)
%       Distances between each wall.
% - path_geometry : struct (no_walls, 1)
%       Set of geometry which the ray will interact with. Does not include
%       the probe or the scatterer.
%
% OUTPUTS:
% - alpha_beta : array (no_legs, 4)
%       The incident and outgoing angles made for each leg. In the 2nd
%       axis, the order of angles provided is:
%           [conv_inc_ang, conv_out_ang, sign_inc_ang, sign_out_ang]
%       where 'conv' is conventional, and 'sign' is signed, both names
%       co-opted from `arim`. In the 1st axis, the angles provided are the
%       incident angles on each wall and the scatterer, and the outgoing
%       angles from the probe and walls.

[no_legs, ~] = size(dists);
alpha_beta = zeros(no_legs, 4);

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
        dist = reshape(-dists(leg, 1:3), 3, 1);
        cart_inc = linsolve(path_geometry(leg).basis, dist);
    else
        cart_inc = -dists(leg, 1:3);
    end
    
    rad_inc = dists(leg, 4);
    % Compute angle.
    alpha = acos(cart_inc(3) / rad_inc);
    % Work out how to shift signed angle.
    azimuth = atan2(cart_inc(2), cart_inc(1));
    if (-pi/2 < azimuth) && (azimuth <= alpha)
        sign_alpha = alpha;
    else
        sign_alpha = -alpha;
    end
    % Work out how to shift conventional ang
    if cart_inc(3) < 0
        conv_alpha = pi - alpha;
    else
        conv_alpha = alpha;
    end
    
    alpha_beta(leg, 1) = conv_alpha;
    alpha_beta(leg, 3) = sign_alpha;
    
    %% Outgoing angles
    
    % Get the vector representing this leg in the wall's coordinate system.
    if leg ~= 1
        dist = reshape(-dists(leg, 1:3), 3, 1);
        cart_out = linsolve(path_geometry(leg-1).basis, dist);
    else
        cart_out = -dists(leg, 1:3);
    end
    
    rad_out = dists(leg, 4);
    % Compute angle.
    beta = acos(cart_out(3) / rad_out);
    % Work out how to shift signed angle.
    azimuth = atan2(cart_out(2), cart_out(1));
    if (-pi/2 < azimuth) && (azimuth <= beta)
        sign_beta = beta;
    else
        sign_beta = -beta;
    end
    % Work out how to shift conventional ang
    if cart_inc(3) < 0
        conv_beta = pi - beta;
    else
        conv_beta = beta;
    end
    
    alpha_beta(leg, 2) = conv_beta;
    alpha_beta(leg, 4) = sign_beta;
    
end

end
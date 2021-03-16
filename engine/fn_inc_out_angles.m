function alpha_beta = fn_inc_out_angles(dists, wall_ids)
% Computes the incident angles (alpha) and outgoing angles (beta) on each
% wall. Assumes that the loop for probes and scatterers is in the parent
% function.

% INPUTS:
% - dists : array of size (no_legs, 4)
%       Distances between each wall.
% - wall_ids : array of size (no_walls+2, 1)
%       Index conveying information on whether each wall the ray interfaces
%       with is a frontwall, backwall, sidewall, scatterer or probe. This 
%       is used to determine how to compute the incident and outgoing 
%       angles.
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
    % If current wall is a frontwall, backwall or scatterer.
    if or(wall_ids(leg+1) == 1, or(wall_ids(leg+1) == 2, ...
          or(wall_ids(leg+1) == -2, wall_ids(leg+1) == -1))) 
        
        % In arim, angles are worked out from the coordinates of the source
        % point of the ray (or the end of the ray in the outgoing case).
        % Here, dists is the direction of the ray. To get the source
        % points, do source = current - ray_direction. To get the angle, we
        % can treat `current` as the origin, so simply take the -ve of the
        % ray direction.
        cart_inc = -dists(leg, 1:3);
        rad_inc = dists(leg, 4);
        alpha = acos(cart_inc(3) / rad_inc);
        azimuth = atan2(cart_inc(2), cart_inc(1));
        if (-pi/2 < azimuth) && (azimuth <= alpha)
            sign_alpha = alpha;
        else
            sign_alpha = -alpha;
        end
        
%         conv_alpha = alpha;
        if cart_inc(3) < 0
            conv_alpha = pi - alpha;
        else
            conv_alpha = alpha;
        end

    elseif wall_ids(leg+1) == 3 % If wall is a sidewall.
        cart_inc = -dists(leg, 1:3);
        rad_inc = dists(leg, 4);
        alpha = acos(cart_inc(1) / rad_inc);
        azimuth = atan2(cart_inc(2), cart_inc(3));
        if (-pi/2 < azimuth) && (azimuth <= alpha)
            sign_alpha = alpha;
        else
            sign_alpha = -alpha;
        end
        
%         conv_alpha = alpha;
        if cart_inc(1) < 0
            conv_alpha = pi - alpha;
        else
            conv_alpha = alpha;
        end

    else % Wall type is not valid.
        disp('Invalid wall type');
    end
    
    alpha_beta(leg, 1) = conv_alpha;
    alpha_beta(leg, 3) = sign_alpha;
end

% Compute outgoing angles from probe and walls.
for leg = 1:no_legs
    % If current wall is a frontwall, backwall or probe.
    if or(wall_ids(leg) == 1, or(wall_ids(leg) == 2, ...
          or(wall_ids(leg) == -1, wall_ids(leg) == -2)))
        
        % Note NDT defn. of angle is between -ve ray direction and +ve
        % axis.
        cart_out = -dists(leg, 1:3);
        rad_out = dists(leg, 4);
        beta = acos(cart_out(3) / rad_out);
        azimuth = atan2(cart_inc(2), cart_inc(1));
        if (-pi/2 < azimuth) && (azimuth <= beta)
            sign_beta = beta;
        else
            sign_beta = -beta;
        end
        
%         conv_beta = beta;
        if cart_inc(3) < 0
            conv_beta = pi - beta;
        else
            conv_beta = beta;
        end
        
    elseif wall_ids(leg) == 3 % If wall is a sidewall.
        cart_out = -dists(leg, 1:3);
        rad_out = dists(leg, 4);
        beta = acos(cart_out(1) / rad_out);
        azimuth = atan2(cart_inc(2), cart_inc(3));
        if (-pi/2 < azimuth) && (azimuth <= beta)
            sign_beta = beta;
        else
            sign_beta = -beta;
        end
        
%         conv_beta = beta;
        if cart_inc(1) < 0
            conv_beta = pi - beta;
        else
            conv_beta = beta;
        end
        
    else % Wall type is not valid.
        disp('Invalid wall type');
    end
    
    alpha_beta(leg, 2) = conv_beta;
    alpha_beta(leg, 4) = sign_beta;
end

end
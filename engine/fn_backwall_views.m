function Views = fn_backwall_views(probe_coords, probe_angle, front_wall, back_wall, mat_speeds, densities, probe_freq, probe_pitch)
% Computes the views arising from reflections from the back wall of the
% geometry. Note that this is an optional function, as these are not always
% desired.
%
% INPUTS:
% - probe_coords : array (probe_els, 3)
%       3D coordinates of the each probe element in the array.
% - probe_angle : double
%       Angle which the probe makes with respect to the geometry in
%       radians. Note that angle is zero when the probe is parallel with
%       the front wall of the geometry (and must always be zero in the
%       contact case).
% - front_wall : array (front_wall_points, 3)
%       3D coordinates of the front wall. Discretised into a number of
%       points.
% - back_wall : array (back_wall_points, 3)
%       3D coordinates of the discretised points of the back wall.
% - mat_speeds : array (3, 1)
%       Speeds of the ray in each medium and for each mode. The order must
%       be [couplant_speed, solid_long_speed, solid_shear_speed]. Note that
%       all three must be provided, even in the contact case.
% - densities : array (2, 1)
%       Density of each medium which the ray passes through. The order must
%       be [couplant_density, solid_density]. Note that couplant density
%       must be provided, even in the contact case.
% - probe_freq : double
%       The frequency of the probe. Note that this is just the probe
%       frequency, not the frequency array which would be computed in the
%       multi-frequency model.
% - probe_pitch : double
%       The pitch of the probe array (i.e. the length of each probe element
%       in the plane being inspected).
%
% OUTPUTS:
% - Views : struct (4, 1)
%       Structure containing the four views generated from the backwall.
%       The order returned is 'L-L', 'L-T', 'T-L', 'T-T'.

% Set up variables required for calculations.
probe_as_scatterer.image_block = probe_coords;
probe_as_scatterer.x_shape = 1;
probe_as_scatterer.z_shape = 1;
probe_as_scatterer.type = "image";

couplant_spd = mat_speeds(1);
mat_long_spd = mat_speeds(2);
mat_shear_spd = mat_speeds(3);

name = ["L-L", "L-T", "T-L", "T-T"];
mat_1_spd = [mat_long_spd, mat_long_spd, mat_shear_spd, mat_shear_spd];
mat_2_spd = [mat_long_spd, mat_shear_spd, mat_long_spd, mat_shear_spd];
modes = [[0, 0]; [0, 1]; [1, 0]; [1, 1]];

% If we are in the contact case, there is no front wall to consider. This
% impacts on the ray weights calculations.
if front_wall == 0
    
    view_speeds = zeros(4, 2); % Shape (no_views, no_legs)
    view_modes  = zeros(4, 2);
    
    for view = 1:4
        view_speeds(view, 1) = mat_1_spd(view);
        view_speeds(view, 2) = mat_2_spd(view);
        view_modes(view, 1)  = modes(view, 1);
        view_modes(view, 2)  = modes(view, 2);
    end
    
    walls = [2];
    medium_ids = [1, 1];
    
    geometry = zeros(1, size(back_wall, 1), size(back_wall, 2));
    geometry(1, :, :) = back_wall;
    
% We must be in the immersion case. We do need to consider the front wall.
else 
    
    view_speeds = zeros(4, 4); % Shape (no_views, no_legs)
    view_modes  = zeros(4, 4);
    
    for view = 1:4
        view_speeds(view, :) = [couplant_spd, mat_1_spd(view), mat_2_spd(view), couplant_spd];
        view_modes(view, :)  = [1, modes(view, 1), modes(view, 2), 1];
    end
    
    walls = [1, 2, 1];
    medium_ids = [0, 1, 1, 0];
    
    geometry = zeros(3, size(front_wall, 1), size(front_wall, 2));
    geometry(1, :, :) = front_wall;
    geometry(2, :, :) = back_wall;
    geometry(3, :, :) = front_wall;
    
end

% Initialise the views.
Views = repmat(struct('name', 'L-L'), 4, 1);

% Compute the ray for each view.
for view = 1:4
    backwall_path = fn_path_info( ...
        name(view), view_modes(view, :), geometry, view_speeds(view, :), mat_speeds, ...
        walls, medium_ids, densities, probe_freq, probe_pitch, probe_coords ...
    );

    ray = fn_compute_ray(probe_as_scatterer, backwall_path, probe_freq);

    [probe_els, ~] = size(ray.min_times);

    Views(view).name = name(view);
    Views(view).min_times = zeros(probe_els ^ 2, 1);
    Views(view).probe_txrx = zeros(probe_els ^ 2, 2);
    Views(view).ray = ray;

    el = 1;

    for t_el = 1 : probe_els
        for r_el = 1 : probe_els

            Views(view).min_times(el, 1) = ray.min_times(t_el, r_el);
            Views(view).probe_txrx(el, :) = [t_el, r_el];

            el = el+1;
        end
    end

    [~, ~, num_freqs] = size(ray.weights.weights);
    Views(view).weights = zeros(probe_els^2, num_freqs);
    el = 1;
    for ii = 1 : probe_els
        for jj = 1 : probe_els
            Views(view).weights(el, :) = ( ...
                ray.weights.weights(ii, jj, :) * ray.weights.inv_directivity(jj, ii, :) ...
            );
        
        el = el+1;
        end
    end
end

end
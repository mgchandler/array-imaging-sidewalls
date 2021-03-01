function fn_tfm(model_config, model_options)
% Computes the sensitivity maps for an array in contact with the solid
% block being inspected. Currently works for a rectangular block with some
% depth defined as the zmax location, and sidewall location defined as the
% xmax location.
%
% INPUTS:
% - model_config : struct (1, 1)
%       A stucture which details the broad config with which to run the
%       model. The fields are detailed below:
%     - PITCH : double
%           The pitch of the probe array. Must be positive real.
%     - PIXEL : integer
%           The pixel resolution of the resulting image.
%     - WALLS : integer
%           The number of points which the walls will be discretised into.
%     - VIEWS : integer
%           A number which specifies which views to model. 1 will image
%           direct views only (total of three imaged); 2 images up to
%           full-skip views from the backwall (total of 21 imaged); 3
%           images up to full-skip views from the sidewall (total of 21
%           imaged); 4 images up to full-skip views from both the sidewall
%           and backwall (total of 55 imaged).
%     - GEOM : logical
%           Logical switch as to whether signals from the geometry will be
%           modelled. In sensitivity, usually this is switched off (i.e.
%           set to 0).
%     - SETUP : logical
%           Logical switch as to whether contact or immersion setup will be
%           taken (0 is contact, 1 is immersion).
% - model_options : OPTIONAL struct (1, 1)
%       Optional structure which contains other options which may be
%       desired. All or none of the following fields may be modelled: if
%       the field is not provided, the default behaviour will take place.
%     - probe_els : integer : DEFAULT = 32
%           The number of elements in the probe. Note that the probe is
%           centered at the coordinates (0, 0, 0) if there is no standoff.
%     - probe_angle : double : DEFAULT = 0
%           The angle the probe makes with the front wall in radians. Note 
%           that in the contact case, this must be set to zero. In the 
%           mmersion case, the standoff must be at least large enough to
%           account foe the swing angle made by the probe when angle is
%           non-zero.
%     - probe_standoff : double : DEFAULT = 0
%           The standoff of the probe with respect to the front wall. Note
%           that in the contact case, this must be set to zero. In the 
%           mmersion case, the standoff must be at least large enough to
%           account foe the swing angle made by the probe when angle is
%           non-zero.
%     - probe_frequency : double : DEFAULT = 5.0e+6
%           The centre frequency of the probe.
%     - el_length : double : DEFAULT = 0.0
%           The additional length each probe element on top of the probe
%           pitch.
%     - geom_shape : struct (1, 1) : DEFAULT = (-25e-3, 25e-3, 0e-3, 40e-3)
%           Structure containing the rectangular geometry of the solid
%           block to be inspected. Required fields are `xmin`, `xmax`,
%           `zmin`, `zmax`.
%     - multi_freq : logical : DEFAULT = 0
%           Logical switch for whether to use the multi-frequency model or
%           not. Default behaviour is to use the single-frequency model.
%     - material_params : struct : DEFAULT air and alum properties
%           A structure detailing the various properties of the couplant
%           and the inspection block. If included, must contain fields
%           titled `couplant_speed`, `couplant_density`, `solid_long_speed`,
%           `solid_shear_speed`, `solid_density`.
%     - scat_info : struct : DEFAULT sdh located at (0, 22e-3).
%           Details on the scatterer which will be modelled. Note that
%           depending on the type of scatterer, different fields are
%           required. Should be an output from the fn_scat_info function.
%     - boxsize : integer : DEFAULT = 0
%           The total length of the square window drawn around the
%           scatterer when finding the peak amplitude in sensitivity
%           calculations.
%     - savepath : string : DEFAULT = array-imaging-sidewalls/output
%           The path where the .fig and .mat files will be saved to.



tic;

%% ---------------------------------------------------------------------- %
% Unpack model_config and model_options                                   %
% ---------------------------------------------------------------------- %%

PITCH = model_config.PITCH;
PIXEL = model_config.PIXEL;
WALLS = model_config.WALLS;
VIEWS = model_config.VIEWS;
GEOM = model_config.GEOM;
SETUP = model_config.SETUP;

probe_els = model_options.probe_els;
probe_angle = model_options.probe_angle;
probe_standoff = model_options.probe_standoff;
probe_frequency = model_options.probe_frequency;
el_length = model_options.el_length + PITCH;
xmin = model_options.geom_shape.xmin;
xmax = model_options.geom_shape.xmax;
zmin = model_options.geom_shape.zmin;
zmax = model_options.geom_shape.zmax;
couplant_speed = model_options.material_params.couplant_speed;
couplant_density = model_options.material_params.couplant_density;
solid_long_speed = model_options.material_params.solid_long_speed;
solid_shear_speed = model_options.material_params.solid_shear_speed;
solid_density = model_options.material_params.solid_density;
boxsize = PIXEL * model_options.boxsize;
scat_info = model_options.scat_info;
savepath = model_options.savepath;

% Work out the size of the final plot from the VIEW parameter.
if VIEWS == 1
    Number_of_ims = 3;
    plot_x = 3;
    plot_z = 1;
    im_width = 450;
    im_height = 240;
    mode = 'Direct';
elseif or(VIEWS == 2, VIEWS == 3)
    Number_of_ims = 21;
    plot_x = 3;
    plot_z = 7;
    im_width = 450;
    im_height = 1050;
    if VIEWS == 2
        mode = 'Backwall Skip';
    else
        mode = 'Sidewall Skip';
    end
else
    Number_of_ims = 55;
    plot_x = 5;
    plot_z = 11;
    im_width = 580;
    im_height = 1280;
    mode = 'Back-Side Skip';
end

% Additional parameters not directly dependent on inputs.
oversampling = 10;
no_cycles = 5;

xsize = xmax - xmin;
zsize = zmax - zmin;

UC = 1e3; % Unit conversion

% Geometry parameters. Do not overlap image params.
front_wall_pos = zmin - 1.0e-5;
front_wall_pixels = WALLS;
back_wall_pos = zmax + 1.0e-5;
back_wall_pixels = WALLS;
side_wall_pos = xmax + 1.0e-5;
side_wall_pixels = WALLS;



%% ------------------------------------------------------------------------
% Set up the scene.
% -------------------------------------------------------------------------

rot_matrix = [cos(probe_angle) 0 sin(probe_angle); 0 1 0; -sin(probe_angle) 0 cos(probe_angle)];

probe_coords = zeros(3, probe_els);
probe_coords(1, :) = linspace(0, (probe_els - 1) * el_length, probe_els);
probe_coords(1, :) = probe_coords(1, :) - mean(probe_coords(1, :));
probe_coords = probe_coords.' * rot_matrix;
probe_coords(:, 3) = probe_coords(:, 3) - probe_standoff;

front_wall = zeros(front_wall_pixels, 3);
front_wall(:, 1) = linspace(xmin, xmax, front_wall_pixels);
front_wall(:, 3) = linspace(front_wall_pos, front_wall_pos, front_wall_pixels);

back_wall = zeros(back_wall_pixels, 3);
back_wall(:, 1) = linspace(xmin, xmax, front_wall_pixels);
back_wall(:, 3) = linspace(back_wall_pos, back_wall_pos, back_wall_pixels);

side_wall = zeros(side_wall_pixels, 3);
side_wall(:, 1) = linspace(side_wall_pos, side_wall_pos, side_wall_pixels);
side_wall(:, 3) = linspace(zmin, zmax, side_wall_pixels);

clear front_wall_pos front_wall_pixels back_wall_pos back_wall_pixels side_wall_pos side_wall_pixels


        
%% ------------------------------------------------------------------------
% Input signal
% -------------------------------------------------------------------------

% Create input signal.
time_step = 1 / (probe_frequency * oversampling); % What is the meaning of oversampling here?
max_t = 1.1 * 4 * (sqrt(xsize ^ 2 + zsize ^ 2) / min(solid_long_speed, solid_shear_speed) + ...
                  (probe_standoff + el_length) / couplant_speed);
time_pts = ceil(max_t / time_step);
[~, ~, freq, in_freq_spec, fft_pts] = fn_create_input_signal(time_pts, probe_frequency, time_step , no_cycles);

time_1 = double(toc);
fprintf('Setup time %.2f secs\n', time_1);
% fprintf('%.2g MB\n', monitor_memory_whos);

clear oversampling no_cycles time_step max_t time in_time_sig



%% ------------------------------------------------------------------------
% Scatterer Simulation path info
% -------------------------------------------------------------------------

tic;

% If we are in the contact case
if ~SETUP
    geometry = 0;

    L_path_info = fn_path_info( ...
        "L", ...
        [0], ...
        geometry, ...
        [solid_long_speed], ...
        [couplant_speed, solid_long_speed, solid_shear_speed], ...
        0, ...
        [1], ...
        [couplant_density, solid_density], ...
        probe_frequency, ...
        el_length, ...
        probe_coords ...
    );
    T_path_info = fn_path_info( ...
        "T", ...
        [1], ...
        geometry, ...
        [solid_shear_speed], ...
        [couplant_speed, solid_long_speed, solid_shear_speed], ...
        0, ...
        [1], ...
        [couplant_density, solid_density], ...
        probe_frequency, ...
        el_length, ...
        probe_coords ...
    ); %#ok<*NBRAK>

    clear geometry

    % -------------------------------------------------------------------------
    if or(VIEWS == 2, VIEWS == 4)
        geometry_bw_skip = zeros(1, size(front_wall, 1), size(front_wall, 2));
        geometry_bw_skip(1, :, :) = back_wall;

        LBL_path_info = fn_path_info( ...
            "LBL", ...
            [0, 0], ...
            geometry_bw_skip, ...
            [solid_long_speed, solid_long_speed], ...
            [couplant_speed, solid_long_speed, solid_shear_speed], ...
            [2], ...
            [1, 1], ...
            [couplant_density, solid_density], ...
            probe_frequency, ...
            el_length, ...
            probe_coords ...
        );
        LBT_path_info = fn_path_info( ...
            "LBT", ...
            [0, 1], ...
            geometry_bw_skip, ...
            [solid_long_speed, solid_shear_speed], ...
            [couplant_speed, solid_long_speed, solid_shear_speed], ...
            [2], ...
            [1, 1], ...
            [couplant_density, solid_density], ...
            probe_frequency, ...
            el_length, ...
            probe_coords ...
        );
        TBL_path_info = fn_path_info( ...
            "TBL", ...
            [1, 0], ...
            geometry_bw_skip, ...
            [solid_shear_speed, solid_long_speed], ...
            [couplant_speed, solid_long_speed, solid_shear_speed], ...
            [2], ...
            [1, 1], ...
            [couplant_density, solid_density], ...
            probe_frequency, ...
            el_length, ...
            probe_coords ...
        );
        TBT_path_info = fn_path_info( ...
            "TBT", ...
            [1, 1], ...
            geometry_bw_skip, ...
            [solid_shear_speed, solid_shear_speed], ...
            [couplant_speed, solid_long_speed, solid_shear_speed], ...
            [2], ...
            [1, 1], ...
            [couplant_density, solid_density], ...
            probe_frequency, ...
            el_length, ...
            probe_coords ...
        );

        clear geometry_bw_skip
    end

    % -------------------------------------------------------------------------
    if or(VIEWS == 3, VIEWS == 4)
        geometry_sw_skip = zeros(1, size(front_wall, 1), size(front_wall, 2));
        geometry_sw_skip(1, :, :) = side_wall;

        LSL_path_info = fn_path_info( ...
            "LSL", ...
            [0, 0], ...
            geometry_sw_skip, ...
            [solid_long_speed, solid_long_speed], ...
            [couplant_speed, solid_long_speed, solid_shear_speed], ...
            [3], ...
            [1, 1], ...
            [couplant_density, solid_density], ...
            probe_frequency, ...
            el_length, ...
            probe_coords ...
        );
        LST_path_info = fn_path_info( ...
            "LST", ...
            [0, 1], ...
            geometry_sw_skip, ...
            [solid_long_speed, solid_shear_speed], ...
            [couplant_speed, solid_long_speed, solid_shear_speed], ...
            [3], ...
            [1, 1], ...
            [couplant_density, solid_density], ...
            probe_frequency, ...
            el_length, ...
            probe_coords ...
        );
        TSL_path_info = fn_path_info( ...
            "TSL", ...
            [1, 0], ...
            geometry_sw_skip, ...
            [solid_shear_speed, solid_long_speed],  ...
            [couplant_speed, solid_long_speed, solid_shear_speed], ...
            [3], ...
            [1, 1], ...
            [couplant_density, solid_density], ...
            probe_frequency, ...
            el_length, ...
            probe_coords ...
        );
        TST_path_info = fn_path_info( ...
            "TST", ...
            [1, 1], ...
            geometry_sw_skip, ...
            [solid_shear_speed, solid_shear_speed], ...
            [couplant_speed, solid_long_speed, solid_shear_speed], ...
            [3], ...
            [1, 1], ...
            [couplant_density, solid_density], ...
            probe_frequency, ...
            el_length, ...
            probe_coords ...
        );

        clear geometry_sw_skip
    end

% If we are in the immersion case.
elseif SETUP

    geometry = zeros(1, size(front_wall, 1), size(front_wall, 2));
    geometry(:, :, :) = front_wall;

    L_path_info = fn_mc_path_info( ...
        "L", ...
        [0, 0], ...
        geometry, ...
        [couplant_speed, solid_long_speed], ...
        [couplant_speed, solid_long_speed, solid_shear_speed], ...
        [1], ...
        [0, 1, 0], ...
        [couplant_density, solid_density], ...
        probe_frequency, ...
        el_length, ...
        probe_coords ...
    );
    T_path_info = fn_mc_path_info( ...
        "T", ...
        [0, 1], ...
        geometry, ...
        [couplant_speed, solid_shear_speed], ...
        [couplant_speed, solid_long_speed, solid_shear_speed], ...
        [1], ...
        [0, 1, 0], ...
        [couplant_density, solid_density], ...
        probe_frequency, ...
        el_length, ...
        probe_coords ...
    );

    clear geometry

    % -------------------------------------------------------------------------
    if or(VIEWS == 2, VIEWS == 4)
        geometry_bw_skip = zeros(2, size(front_wall, 1), size(front_wall, 2));
        geometry_bw_skip(1, :, :) = front_wall;
        geometry_bw_skip(2, :, :) = back_wall;

        LBL_path_info = fn_mc_path_info( ...
            "LBL", ...
            [0, 0, 0], ...
            geometry_bw_skip, ...
            [couplant_speed, solid_long_speed, solid_long_speed], ...
            [couplant_speed, solid_long_speed, solid_shear_speed], ...
            [1, 2], ...
            [0, 1, 1, 0], ...
            [couplant_density, solid_density], ...
            probe_frequency, ...
            el_length, ...
            probe_coords ...
        );
        LBT_path_info = fn_mc_path_info( ...
            "LBT", ...
            [0, 0, 1], ...
            geometry_bw_skip, ...
            [couplant_speed, solid_long_speed, solid_shear_speed], ...
            [couplant_speed, solid_long_speed, solid_shear_speed], ...
            [1, 2], ...
            [0, 1, 1, 0], ...
            [couplant_density, solid_density], ...
            probe_frequency, ...
            el_length, ...
            probe_coords ...
        );
        TBL_path_info = fn_mc_path_info( ...
            "TBL", ...
            [0, 1, 0], ...
            geometry_bw_skip, ...
            [couplant_speed, solid_shear_speed, solid_long_speed], ...
            [couplant_speed, solid_long_speed, solid_shear_speed], ...
            [1, 2], ...
            [0, 1, 1, 0], ...
            [couplant_density, solid_density], ...
            probe_frequency, ...
            el_length, ...
            probe_coords ...
        );
        TBT_path_info = fn_mc_path_info( ...
            "TBT", ...
            [0, 1, 1], ...
            geometry_bw_skip, ...
            [couplant_speed, solid_shear_speed, solid_shear_speed], ...
            [couplant_speed, solid_long_speed, solid_shear_speed], ...
            [1, 2], ...
            [0, 1, 1, 0], ...
            [couplant_density, solid_density], ...
            probe_frequency, ...
            el_length, ...
            probe_coords ...
        );

        clear geometry_bw_skip
    end

    % -------------------------------------------------------------------------
    if or(VIEWS == 3, VIEWS == 4)
        geometry_sw_skip = zeros(2, size(front_wall, 1), size(front_wall, 2));
        geometry_sw_skip(1, :, :) = front_wall;
        geometry_sw_skip(2, :, :) = side_wall;

        LSL_path_info = fn_mc_path_info( ...
            "LSL", ...
            [0, 0, 0], ...
            geometry_sw_skip, ...
            [couplant_speed, solid_long_speed, solid_long_speed], ...
            [couplant_speed, solid_long_speed, solid_shear_speed], ...
            [1, 3], ...
            [0, 1, 1, 0], ...
            [couplant_density, solid_density], ...
            probe_frequency, ...
            el_length, ...
            probe_coords ...
        );
        LST_path_info = fn_mc_path_info( ...
            "LST", ...
            [0 0 1], ...
            geometry_sw_skip, ...
            [couplant_speed, solid_long_speed, solid_shear_speed], ...
            [couplant_speed, solid_long_speed, solid_shear_speed], ...
            [1 3], ...
            [0, 1, 1, 0], ...
            [couplant_density, solid_density], ...
            probe_frequency, ...
            el_length, ...
            probe_coords ...
        );
        TSL_path_info = fn_mc_path_info( ...
            "TSL", ...
            [0, 1, 0], ...
            geometry_sw_skip, ...
            [couplant_speed, solid_shear_speed, solid_long_speed], ...
            [couplant_speed, solid_long_speed, solid_shear_speed], ...
            [1, 3], ...
            [0, 1, 1, 0], ...
            [couplant_density, solid_density], ...
            probe_frequency, ...
            el_length, ...
            probe_coords ...
        );
        TST_path_info = fn_mc_path_info( ...
            "TST", ...
            [0, 1, 1], ...
            geometry_sw_skip, ...
            [couplant_speed, solid_shear_speed, solid_shear_speed], ...
            [couplant_speed, solid_long_speed, solid_shear_speed], ...
            [1, 3], ...
            [0, 1, 1, 0], ...
            [couplant_density, solid_density], ...
            probe_frequency, ...
            el_length, ...
            probe_coords ...
        );

        clear geometry_sw_skip
    end
    
end

time_2 = double(toc);

fprintf('Scatterer paths computed in %.2f secs\n', time_2);
% fprintf('%.2g MB\n', monitor_memory_whos);

clear couplant_speed couplant_density mat_long_speed mat_shear_speed mat_density



%% ------------------------------------------------------------------------
% Views
% -------------------------------------------------------------------------

xpts = round(xsize / PIXEL);
zpts = round(zsize / PIXEL);

% Compute Imaging Paths and Views
L_path_scat = fn_compute_ray(scat_info, L_path_info, probe_frequency);
T_path_scat = fn_compute_ray(scat_info, T_path_info, probe_frequency);
Paths = [L_path_scat T_path_scat];
Names = ["L" "T"];
Names_rev = ["L" "T"];
clear L_path_scat T_path_scat X Z
if or(VIEWS == 2, VIEWS == 4)
    LBL_path_scat = fn_compute_ray(scat_info, LBL_path_info, probe_frequency);
    LBT_path_scat = fn_compute_ray(scat_info, LBT_path_info, probe_frequency);
    TBL_path_scat = fn_compute_ray(scat_info, TBL_path_info, probe_frequency);
    TBT_path_scat = fn_compute_ray(scat_info, TBT_path_info, probe_frequency);
    Paths = [Paths LBL_path_scat LBT_path_scat TBL_path_scat TBT_path_scat];
    Names = [Names "LBL" "LBT" "TBL" "TBT"];
    Names_rev = [Names_rev "LBL" "TBL" "LBT" "TBT"];
    clear LBL_path_scat LBT_path_scat TBL_path_scat TBT_path_scat
end
if or(VIEWS == 3, VIEWS == 4)
    LSL_path_scat = fn_compute_ray(scat_info, LSL_path_info, probe_frequency);
    LST_path_scat = fn_compute_ray(scat_info, LST_path_info, probe_frequency);
    TSL_path_scat = fn_compute_ray(scat_info, TSL_path_info, probe_frequency);
    TST_path_scat = fn_compute_ray(scat_info, TST_path_info, probe_frequency);
    Paths = [Paths LSL_path_scat LST_path_scat TSL_path_scat TST_path_scat];
    Names = [Names "LSL" "LST" "TSL" "TST"];
    Names_rev = [Names_rev "LSL" "TSL" "LST" "TST"];
    clear LSL_path_im LST_path_im TSL_path_im TST_path_im
end

Views = repmat(fn_create_view(Paths(1), Paths(1), scat_info), Number_of_ims, 1);
Ims = repmat(fn_create_im("-", xpts+1, zpts+1), Number_of_ims, 1);
Namelist = repmat("-", Number_of_ims, 1);
i = 1;
for t_path = 1 : size(Paths, 2)
    for r_path = 1 : size(Paths, 2)
        name = strcat(Names(t_path), "-", Names_rev(r_path));
        revname = strcat(Names(r_path), "-", Names_rev(t_path));
        if ~any(strcmp(Namelist, revname))
            Namelist(i) = name;
            Views(i) = fn_create_view(Paths(t_path), Paths(r_path), scat_info);
            Views(i).name = name;
            Sens(i).name = name;
            i = i + 1;
        end
    end
end

clear boxsize Namelist



%% ------------------------------------------------------------------------
% Geometry simulation
% -------------------------------------------------------------------------
tic;
if GEOM
    backwall_views = fn_backwall_views( ...
        probe_coords, probe_angle, 0, back_wall, [couplant_speed solid_long_speed solid_shear_speed], ...
        [couplant_density solid_density], probe_frequency, el_length ...
    );
    amp_b = repmat(zeros(1, size(backwall_views(1).min_times, 1)), size(backwall_views, 2), 1);
    for view = 1 : size(backwall_views, 2)
        amp_b(view, :) = conj( ...
            backwall_views(view).directivity1 .* backwall_views(view).directivity2 .* ...
            backwall_views(view).bs1 .* backwall_views(view).tr1 ...
        );
    end
    
    Geo = 'Geo';
else
    Geo = 'NGeo';
end

clear probe_pitch

time_3 = double(toc);

fprintf('Rays traced and views assembled in %.2f secs\n', time_3);
% fprintf('%.2g MB\n', monitor_memory_whos);



%% ------------------------------------------------------------------------
% Simulation
% -------------------------------------------------------------------------

% Scatterer simulation.
tic;
out_freq_spec = 0;

for view = 1 : size(Views, 1)
    ray_weights_path1 = fn_compute_ray_weights(Views(view).path_1, probe_frequency);
    ray_weights_path2 = fn_compute_ray_weights(Views(view).path_2, probe_frequency);
    ray_weights = zeros(probe_els^2, size(scat_info.image_block, 1), 1);
    el = 1;
    for tx = 1 : probe_els
        for rx = 1 : probe_els
            ray_weights(el, :, :) = ray_weights_path1.weights(tx, :, :) .* ray_weights_path2.inv_weights(rx, :, :);
            el = el + 1;
        end
    end
    scat_amp = fn_scattering_amps(Views(view), probe_frequency);
    for scatterer = 1 : size(scat_info.image_block, 1)
        amp = conj(scat_amp(:, scatterer, :) .* ray_weights(:, scatterer, :));
        out_freq_spec = out_freq_spec + fn_propagate_spectrum_mc_2(freq, in_freq_spec, Views(view).min_times(:, scatterer), amp, 0);
    end
end

% Convert back to time.
[FMC.time, FMC.time_data] = fn_convert_spectrum_to_time(freq, out_freq_spec, fft_pts, time_pts);
FMC.time = FMC.time';

% Hilbert Filtering (?) from fast_DAS function
diagonals = spdiags([1:length(FMC.time)]' < length(FMC.time)/2, 0, length(FMC.time), length(FMC.time));
FMC.time_data = ifft(diagonals * fft(FMC.time_data));



clear time_pts freq in_freq_spec fft_pts box_pts xpts_im zpts_im grid_pt scat_pt sens_i sens_k
clear weights scat_amp amp out_freq_spec diagonals tau Im scatterer_coords
% clear image_block_info
if GEOM
    clear backwall_views amp_b
end

time_4 = double(toc);

fprintf('Scatterer simulated in %.2f secs\n', time_4);
% fprintf('%.2g MB\n', monitor_memory_whos);



%% ------------------------------------------------------------------------
% Imaging
% -------------------------------------------------------------------------

tic

db_range_for_output = 40;
xpts = round(xsize / PIXEL);
zpts = round(zsize / PIXEL);
im_x = linspace(xmin, xmax, xpts+1);
im_z = linspace(zmin, zmax, zpts+1);
[X, Z] = meshgrid(im_x, im_z);

pt = 1;
image_block = zeros((xpts+1)*(zpts+1), 3);

for xpt = 1 : xpts+1
    for zpt = 1 : zpts+1
        image_block(pt, 1) = X(zpt, xpt);
        image_block(pt, 3) = Z(zpt, xpt);
        pt = pt + 1;
    end
end

image_block_info = fn_scat_info("image", image_block);



tic;
% Compute Imaging Paths and Views
disp('Computing imaging')

L_path_im = fn_compute_ray(image_block_info, L_path_info);
T_path_im = fn_compute_ray(image_block_info, T_path_info);
Paths_im = [L_path_im, T_path_im];

if or(VIEWS == 2, VIEWS == 4)
    LBL_path_im = fn_compute_ray(image_block_info, LBL_path_info);
    LBT_path_im = fn_compute_ray(image_block_info, LBT_path_info);
    TBL_path_im = fn_compute_ray(image_block_info, TBL_path_info);
    TBT_path_im = fn_compute_ray(image_block_info, TBT_path_info);
    Paths_im = [Paths_im, LBL_path_im, LBT_path_im, TBL_path_im, TBT_path_im];
end
if or(VIEWS == 3, VIEWS == 4)
    LSL_path_im = fn_compute_ray(image_block_info, LSL_path_info);
    LST_path_im = fn_compute_ray(image_block_info, LST_path_info);
    TSL_path_im = fn_compute_ray(image_block_info, TSL_path_info);
    TST_path_im = fn_compute_ray(image_block_info, TST_path_info);
    Paths_im = [Paths_im, LSL_path_im, LST_path_im, TSL_path_im, TST_path_im];
end

tic;
Views_im = repmat(fn_create_view(Paths_im(1), Paths_im(1)), Number_of_ims, 1);
Ims = repmat(fn_create_im("-", xpts+1, zpts+1), Number_of_ims, 1);
Namelist = [];
i = 1;
for t_path = 1 : size(Paths_im, 2)
    for r_path = 1 : size(Paths_im, 2)
        name = strcat(Names(t_path), "-", Names_rev(r_path));
        revname = strcat(Names(r_path), "-", Names_rev(t_path));
        if ~any(strcmp(Namelist, revname))
            Namelist = [Namelist name];
            Views_im(i) = fn_create_view(Paths_im(t_path), Paths_im(r_path));
            Ims(i).name = name;
            i = i + 1;
        end
    end
end

tic;
% Lookup times in FMC data.
for tr_pair = 1 : probe_els ^ 2
    for view = 1 : size(Views_im, 1)
        tau = reshape(Views_im(view).min_times(tr_pair, :), [zpts+1, xpts+1]);
        Ims(view).image = Ims(view).image + interp1(FMC.time, FMC.time_data(:, tr_pair), tau, 'linear', 0);
    end
end

time_5 = double(toc);

fprintf('Scatterer imaged in %.2f secs\n', time_4);
% fprintf('%.2g MB\n', monitor_memory_whos);



tic

max_ = 0;
for view = 1 : size(Views_im, 1)
    if max(abs(Ims(view).image(:))) > max_
        max_ = max(abs(Ims(view).image(:)));
    end
end

for view = 1 : size(Views_im, 1)
   Ims(view).db_image = 20 * log10(abs(Ims(view).image) ./ max_); 
end

if VIEWS == 1
    View_names = ["L-L", "L-T", "T-T"];
elseif VIEWS == 2
    View_names = ["L-L", "L-T", "T-T", "L-LBL", "T-LBL", "L-TBL", "T-TBL", ...
                  "L-LBT", "T-LBT", "L-TBT", "T-TBT", "LBL-LBL", "LBL-LBT", "LBL-TBL", ...
                  "LBL-TBT", "LBT-LBT", "LBT-TBL", "LBT-TBT", "TBL-LBT", "TBL-TBT", "TBT-TBT"];
elseif VIEWS == 3
    View_names = ["L-L", "L-T", "T-T", "L-LSL", "T-LSL", "L-TSL", "T-TSL", ...
                  "L-LST", "T-LST", "L-TST", "T-TST", "LSL-LSL", "LSL-LST", "LSL-TSL", ...
                  "LSL-TST", "LST-LST", "LST-TSL", "LST-TST", "TSL-LST", "TSL-TST", "TST-TST"];
elseif VIEWS == 4
    View_names = ["L-L", "L-T", "T-T", "LBL-T", "LBL-T", "LBT-L", "LBT-T", ...
                  "TBL-L", "TBL-T", "TBT-L", "TBT-T", "LBL-LBL", "LBL-LBT", "LBL-TBL", ...
                  "LBL-TBT", "LBT-LBT", "LBT-TBL", "LBT-TBT", "TBL-LBT", "TBL-TBT", "TBT-TBT", ...
                  "L-LSL", "T-LSL", "L-TSL", "T-TSL", ...
                  "L-LST", "T-LST", "L-TST", "T-TST", "LSL-LSL", "LSL-LST", "LSL-TSL", ...
                  "LSL-TST", "LST-LST", "LST-TSL", "LST-TST", "TSL-LST", "TSL-TST", "TST-TST"];
end

sort_idx = zeros(Number_of_ims, 1);
view_names = struct2cell(Views);
view_names = view_names(3, :);
for view = 1:Number_of_ims
    sort_idx(view) = find(strcmp([view_names{:}], View_names(view)));
end

% Plot.
im_x = linspace(xmin, xmax, xpts+1);
im_z = linspace(zmin, zmax, zpts+1);
fig = figure(1);
ax = repmat(subplot(plot_z, plot_x, 1), Number_of_ims, 1);
for im = 1:Number_of_ims
    ax(im) = subplot(plot_z, plot_x, im);
    imagesc(im_x*UC, im_z*UC, Ims(sort_idx(im)).db_image);
    hold on
    title(Ims(sort_idx(im)).name)
    caxis([-db_range_for_output, 0])
    plot(probe_coords(:, 1)*UC, probe_coords(:, 3)*UC, 'go');
    plot(front_wall(:, 1)*UC, front_wall(:, 3)*UC, 'r');
    plot(back_wall(:, 1)*UC, back_wall(:, 3)*UC, 'r');
    plot(side_wall(:, 1)*UC, side_wall(:, 3)*UC, 'r');
    for s = 1 : size(scat_info.image_block, 1)
        rectangle('Position', [scat_info.image_block(s, 1)*UC - 2.5, scat_info.image_block(s, 3)*UC - 2.5, 5, 5], 'EdgeColor', 'r');
    end
    axis equal; axis tight;
end

% Resize image and move for colorbar
set(fig, 'Position', [20, 20, im_width, im_height])
posa = cell2mat(get(ax, 'Position'));
h = colorbar;
set(ax(im), 'Position', posa(im, :))
set(ax, 'units', 'pix')
set(h, 'units', 'pix')
posf = get(fig, 'Position'); % gives x left, y bottom, width, height
set(fig, 'Position',  [posf(1:2) posf(3)*1.1 posf(4)])
hpos = h.Position;
posa = cell2mat(get(ax, 'Position'));
pos1 = posa(1, :);
set(h, 'Position', [hpos(1)+10, hpos(2), hpos(3)*2, pos1(2)-hpos(2)+pos1(4)])

h.Label.String = 'dB';



cd(savepath)
filename_fig = sprintf('Sens Contact %s %s %s scatterer - p=%.2e pix=%.2e walls=%d - shape=(%.2e,%.2e,%.2e,%.2e).fig', mode, Geo, image_block_info.type, PITCH, PIXEL, WALLS, xmin, xmax, zmin, zmax);
filename_mat = sprintf('Sens Contact %s %s %s scatterer - p=%.2e pix=%.2e walls=%d - shape=(%.2e,%.2e,%.2e,%.2e).mat', mode, Geo, image_block_info.type, PITCH, PIXEL, WALLS, xmin, xmax, zmin, zmax);
saveas(fig, filename_fig)
close all

time_6 = double(toc);

fprintf('Plotted in %.2f secs\n', time_6);
% fprintf('%.2g MB\n', monitor_memory_whos);

times = [time_1, time_2, time_3, time_4, time_5, time_6];

save(filename_mat, 'times', 'Ims', 'Views')

clear

end
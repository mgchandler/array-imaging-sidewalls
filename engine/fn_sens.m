function fn_sens(model_config, model_options)
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
%           modelled. In sensitivity,  this should be switched off (i.e.
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
%     - savename : string : DEFAULT = "sens MODE GEOM VIEWS PITCH PIXEL WALLS"
%           The name which will files will be saved as. The file types
%           (.fig and .mat) will be appended to this string.



tic;

%% ---------------------------------------------------------------------- %
% Unpack model_config and model_options                                   %
% ---------------------------------------------------------------------- %%

PITCH = model_config.PITCH;
PIXEL = model_config.PIXEL;
WALLS = model_config.WALLS;
VIEWS = model_config.VIEWS;
assert(~model_config.GEOM, 'Error in fn_sens model_config.GEOM must be 0 for sensitivity.');
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
boxsize = model_options.boxsize;
image_block_info = model_options.scat_info;
savepath = model_options.savepath;
savename = model_options.savename;
geometry = model_options.geometry;
max_no_reflections = model_options.max_no_reflections;

% Work out the size of the final plot from the VIEW parameter.

% Additional parameters not directly dependent on inputs.
oversampling = 10;
no_cycles = 5;

no_walls = size(geometry, 1);

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

clear PITCH WALLS



%% ---------------------------------------------------------------------- %
% Set up the scene.                                                       %
% ---------------------------------------------------------------------- %%

rot_matrix = [cos(probe_angle) 0 sin(probe_angle); 0 1 0; -sin(probe_angle) 0 cos(probe_angle)];

probe_coords = zeros(3, probe_els);
probe_coords(1, :) = linspace(0, (probe_els - 1) * el_length, probe_els);
probe_coords(1, :) = probe_coords(1, :) - mean(probe_coords(1, :));
probe_coords = probe_coords.' * rot_matrix;
probe_coords(:, 3) = probe_coords(:, 3) - probe_standoff - 1e-5;

clear probe_angle front_wall_pos back_wall_pos side_wall_pos front_wall_pixels back_wall_pixels side_wall_pixels
clear rot_matrix



%% ---------------------------------------------------------------------- %
% Input signal                                                            %
% ---------------------------------------------------------------------- %%

% Create input signal.
time_step = 1 / (probe_frequency * oversampling); % What is the meaning of oversampling here?
max_t = 1.1 * 4 * (sqrt(xsize ^ 2 + zsize ^ 2) / min(solid_long_speed, solid_shear_speed) + ...
                  (probe_standoff + el_length*probe_els) / couplant_speed);
time_pts = ceil(max_t / time_step);
[~, ~, freq, in_freq_spec, fft_pts] = fn_create_input_signal(time_pts, probe_frequency, time_step , no_cycles);

clear probe_standoff oversampling no_cycles time_step max_t



%% ---------------------------------------------------------------------- %
% Scatterer Simulation path info - Contact                                %
% ---------------------------------------------------------------------- %%

tic;
    
mode_names = ["L", "T"];
speeds = [solid_long_speed, solid_shear_speed];

Path_info_list = []; % Need to preallocate this field. Will do later.

% If we are in the contact case
if ~SETUP
    path_geometry = 0;
    
% ----------------------------------------------------------------------- %
% Direct Paths                                                            %
% ----------------------------------------------------------------------- %

    for mode = 0:1
        mode_name = mode_names(mode+1);
        Direct_path_info = fn_path_info( ...
            mode_name, ...
            mode_name, ...
            [mode], ...
            path_geometry, ...
            speeds(mode+1), ...
            [couplant_speed, solid_long_speed, solid_shear_speed], ...
            0, ... % index for wall identity
            [1], ... % index for material identity
            [couplant_density, solid_density], ...
            probe_frequency, ...
            el_length, ...
            probe_coords ...
        );
        Path_info_list = [Path_info_list, Direct_path_info];
        clear Direct_path_info
    end

    clear path_geometry
    
% ----------------------------------------------------------------------- %
% Skip Paths                                                              %
% ----------------------------------------------------------------------- %

    if max_no_reflections > 0
        for wall = 1:no_walls
            path_geometry = geometry(wall);
            for mode1 = 0:1 % Mode of the first leg
                mode1_name = mode_names(mode1+1);
                for mode2 = 0:1 % Mode of the second leg
                    mode2_name = mode_names(mode2+1);
                    Skip_path_info = fn_path_info( ...
                        sprintf("%s %s %s", mode1_name, path_geometry.name, mode2_name), ...
                        sprintf("%s %s %s", mode2_name, path_geometry.name, mode1_name), ...
                        [mode1, mode2], ...
                        path_geometry, ...
                        [speeds(mode1+1), speeds(mode2+1)], ...
                        [couplant_speed, solid_long_speed, solid_shear_speed], ...
                        [path_geometry.wall_id], ...
                        [1, 1], ...
                        [couplant_density, solid_density], ...
                        probe_frequency, ...
                        el_length, ...
                        probe_coords ...
                    );
                    Path_info_list = [Path_info_list, Skip_path_info];
                    clear Skip_path_info
                end
            end
        end
    end

%% ---------------------------------------------------------------------- %
% Scatterer Simulation path info - Contact                                %
% ---------------------------------------------------------------------- %%

% If we are in the immersion case.
elseif SETUP
    
    path_geometry = geometry(1);
    
% ----------------------------------------------------------------------- %
% Direct Paths                                                            %
% ----------------------------------------------------------------------- %

    for mode = 0:1
        mode_name = mode_names(mode+1);
        Direct_path_info = fn_path_info( ...
            mode_name, ...
            mode_name, ...
            [0, mode], ...
            path_geometry, ...
            [couplant_speed, speeds(mode+1)], ...
            [couplant_speed, solid_long_speed, solid_shear_speed], ...
            [path_geometry.wall_id], ...
            [0, 1], ...
            [couplant_density, solid_density], ...
            probe_frequency, ...
            el_length, ...
            probe_coords ...
        );
        Path_info_list = [Path_info_list, Direct_path_info];
        clear Direct_path_info
    end
    
% ----------------------------------------------------------------------- %
% Skip Paths                                                              %
% ----------------------------------------------------------------------- %

    if max_no_reflections > 0
        for wall = 2:no_walls
            path_geometry = repmat(geometry(1), 2, 1);
            path_geometry(2) = geometry(wall);
            for mode1 = 0:1 % Mode of the first leg
                mode1_name = mode_names(mode1+1);
                for mode2 = 0:1 % Mode of the second leg
                    mode2_name = mode_names(mode2+1);
                    Skip_path_info = fn_path_info( ...
                        sprintf("%s %s %s", mode1_name, path_geometry(2).name, mode2_name), ...
                        sprintf("%s %s %s", mode2_name, path_geometry(2).name, mode1_name), ...
                        [0, mode1, mode2], ...
                        path_geometry, ...
                        [couplant_speed, speeds(mode1+1), speeds(mode2+1)], ...
                        [couplant_speed, solid_long_speed, solid_shear_speed], ...
                        [path_geometry(1).wall_id, path_geometry(2).wall_id], ...
                        [0, 1, 1], ...
                        [couplant_density, solid_density], ...
                        probe_frequency, ...
                        el_length, ...
                        probe_coords ...
                    );
                    Path_info_list = [Path_info_list, Skip_path_info];
                    clear Skip_path_info
                end
            end
        end
    end
    
else
    error('Error in fn_sens: Invalid SETUP - must be contact (0) or immersion (1).')
end

time_1 = double(toc);

fprintf('Setup time %.2f secs.\n', time_1);
fprintf('%.2g MB used.\n', monitor_memory_whos);

clear SETUP el_length couplant_speed couplant_density solid_long_speed solid_shear_speed solid_density



%% ---------------------------------------------------------------------- %
% Scatterer and Imaging setup                                             %
% ---------------------------------------------------------------------- %%
% Precompute views for all scatterers, as well as imaging views.

tic;

db_range_for_output = 40;

% The total number of points where scatterers will be placed in the x- and
% z- axes.
xpts = round(xsize / PIXEL) + 2*boxsize;
zpts = round(zsize / PIXEL) + 2*boxsize;
% The number of points in the box plotted for each individual
% scatterer, from which the maximum is taken for the sensitivity plot.
boxpix = boxsize * PIXEL;
% Helper arrays for use in the meshgrid.
x = linspace(xmin-boxpix, xmax+boxpix, xpts+1);
z = linspace(zmin-boxpix, zmax+boxpix, zpts+1);
% Meshgrid for assembling the image block.
[X, Z] = meshgrid(x, z);
% [im_X, im_Z] = meshgrid(im_x, im_z);

pt = 0;
image_block = zeros((xpts+1)*(zpts+1), 3);

for xpt = 1 : xpts+1
    for zpt = 1 : zpts+1
        pt = pt + 1;
        image_block(pt, 1) = X(zpt, xpt);
        image_block(pt, 3) = Z(zpt, xpt);
    end
end

% Reuse imaging paths for 
image_block_info.image_block = image_block;
scatterer_coords = reshape(image_block, zpts+1, xpts+1, 3);

% Compute Imaging Paths
Paths = [];
for path = 1:size(Path_info_list, 2)
    Paths = [Paths, fn_compute_ray(image_block_info, Path_info_list(path), probe_frequency)];
end
clear Path_info_list X Z pt image_block

% Create views from these paths.
Views = fn_make_views(VIEWS, Paths, 1);
clear Paths
Number_of_ims = size(Views, 1);
if Number_of_ims == 3
    plot_x = 3;
    plot_z = 1;
    im_width = 450;
    im_height = 240;
elseif Number_of_ims == 21
    plot_x = 3;
    plot_z = 7;
    im_width = 450;
    im_height = 1050;
elseif Number_of_ims == 55
    plot_x = 5;
    plot_z = 11;
    im_width = 580;
    im_height = 1280;
else
    err_msg = sprintf('fn_sens: Unexpected number of images being plotted.\n%d image(s) being plotted.', Number_of_ims);
    error(err_msg)
end
Sens = repmat(fn_create_im("-", xpts+1, zpts+1), Number_of_ims, 1);
for view = 1 : Number_of_ims
    Sens(view).name = Views(view).name;
end

time_2 = double(toc);

fprintf('Rays traced in %.2f secs.\n', time_2);
fprintf('%.2g MB used.\n', monitor_memory_whos);

clear probe_frequency boxpix



%% ---------------------------------------------------------------------- %
% Sensitivity Images                                                      %
% ---------------------------------------------------------------------- %%

tic

% Obscure points outside of the geometry: these will be zero in logical
% checks below, so multiply out all checks to kill pixels introduced by
% boxsize variable.
are_points_in_geometry = (scatterer_coords(:,:,1) >= xmin) .* ...
                         (scatterer_coords(:,:,1) <= xmax) .* ...
                         (scatterer_coords(:,:,3) >= zmin) .* ...
                         (scatterer_coords(:,:,3) < zmax);
                     
% Set up sens_i and sens_k indices for referencing when setting values of
% the Sens image.
sens_i = [1:xpts+1]+xpts+1;
sens_k = [1:zpts+1]+zpts+1;
sens_i_min = sens_i-boxsize;
sens_i_max = sens_i+boxsize;
sens_k_min = sens_k-boxsize;
sens_k_max = sens_k+boxsize;

clear xmin xmax zmin zmax sens_i sens_k scatterer_coords



%% ---------------------------------------------------------------------- %
% Start of sensitivity loop                                               %
% ---------------------------------------------------------------------- %%

tic;

grid_pt = 0;
for xpt_im = 1:xpts+1
    for zpt_im = 1:zpts+1
        grid_pt = grid_pt + 1;
        
        

% ----------------------------------------------------------------------- %
% Simulation Step                                                         %
% ----------------------------------------------------------------------- %

        for view = 1 : Number_of_ims
            weights = Views(view).weights(:, grid_pt, 1);
            scat_amp = Views(view).scat_amps(:, grid_pt, 1);
            amp = conj(scat_amp .* weights);
            out_freq_spec = fn_propagate_spectrum_mc(freq, in_freq_spec, Views(view).min_times(:, grid_pt), amp, 0);
            
            clear weights scat_amp amp

            % Convert back to time.
            [FMC_time, FMC_time_data] = fn_convert_spectrum_to_time(freq, out_freq_spec, fft_pts, time_pts);
            FMC_time = FMC_time';
            
            clear out_freq_spec

            % Hilbert Filtering (?) from fast_DAS function
            diagonals = spdiags([1:length(FMC_time)]' < length(FMC_time)/2, 0, length(FMC_time), length(FMC_time));
            FMC_time_data = ifft(diagonals * fft(FMC_time_data));
            
            clear diagonals
        
% ----------------------------------------------------------------------- %
% Imaging Step                                                            %
% ----------------------------------------------------------------------- %
            
            if boxsize == 0
                tau = reshape(Views(view).min_times, [probe_els^2, zpts+1, xpts+1]);
                Im = (are_points_in_geometry(zpt_im, xpt_im) * ...
                    sum(diag(interp1(FMC_time, FMC_time_data, tau(:, zpt_im, xpt_im), 'linear', 0))) ...
                );
            else
                tau = repmat(reshape(Views(view).min_times, [probe_els^2, zpts+1, xpts+1]),1,3,3);
                Im = zeros(boxsize*2+1);
                for s_i = sens_i_min(xpt_im):sens_i_max(xpt_im)
                    for s_k = sens_k_min(zpt_im):sens_k_max(zpt_im)
                        Im(s_k-sens_k_min(zpt_im)+1, s_i-sens_i_min(xpt_im)+1) = ( ...
                            are_points_in_geometry(zpt_im, xpt_im) * ...
                            sum(diag(interp1(FMC_time, FMC_time_data, tau(:, s_k, s_i), 'linear', 0))) ...
                        );
                    end
                end
            end
%             tau2 = tau(:, sens_k_min(zpt):sens_k_max(zpt), sens_i_min(xpt):sens_i_max(xpt));
%             interp = interp1(FMC_time, FMC_time_data, tau(:, sens_k_min(zpt):sens_k_max(zpt), sens_i_min(xpt):sens_i_max(xpt)), 'linear', 0);
%             Im = sum(diag(interp1(FMC_time, FMC_time_data, tau(:, sens_k_min(zpt):sens_k_max(zpt), sens_i_min(xpt):sens_i_max(xpt)), 'linear', 0)));
            
            Sens(view).image(zpt_im, xpt_im) = max(Im, [], 'all');
            
            clear FMC_time FMC_time_data tau Im

        end
        
%% ---------------------------------------------------------------------- %
% Estimate Loop Runtime                                                   %
% ---------------------------------------------------------------------- %%

        if and(xpt_im == 1, zpt_im == 1)
            timing = double(toc);
            fprintf('    One loop performed in %.1f secs\n', timing);
            est_runtime = timing * xsize * zsize / (PIXEL^2);
            if est_runtime < 60
                fprintf('    Estimated runtime is %.1f secs\n', est_runtime);
            elseif (est_runtime / 60) < 60
                fprintf('    Estimated runtime is %.1f mins\n', est_runtime/60);
            else
                fprintf('    Estimated runtime is %.2f hrs\n', est_runtime/3600);
            end
        elseif and(xpt_im == round((xpts+1)/2), zpt_im == 1)
            timing = double(toc);
            est_runtime = timing * 2;
            if est_runtime < 60
                fprintf('    Updated est runtime is %.1f secs\n', est_runtime);
            elseif (est_runtime / 60) < 60
                fprintf('    Updated est runtime is %.1f mins\n', est_runtime/60);
            else
                fprintf('    Updated est runtime is %.2f hrs\n', est_runtime/3600);
            end
        end

    end
end

time_3 = double(toc);

fprintf('Finished Sensitivity Loop in %.2f secs\n', time_3);
fprintf('%.2g MB used.\n', monitor_memory_whos);

clear PIXEL probe_els boxsize xsize zsize time_pts freq in_freq_spec fft_pts
clear are_points_in_geometry sens_i_min sens_i_max sens_k_min sens_k_max grid_pt



%% ---------------------------------------------------------------------- %
% Convert to dB and Plot                                                  %
% ---------------------------------------------------------------------- %%
tic;

max_ = 0;
for view = 1 : size(Views, 1)
    if max(abs(Sens(view).image(:))) > max_
        max_ = max(abs(Sens(view).image(:)));
    end
end

for view = 1 : size(Views, 1)
   Sens(view).db_image = 20 * log10(abs(Sens(view).image) ./ max_); 
end

% Plot.
fig = figure(1);
ax = repmat(subplot(plot_z, plot_x, 1), Number_of_ims, 1);
if image_block_info.type == "crack"
    sgtitle(sprintf('Sens %.2f Crack - %.2f deg', image_block_info.crack_length, rad2deg(image_block_info.angle)))
end
for im = 1:Number_of_ims
    ax(im) = subplot(plot_z, plot_x, im);
    imagesc(x*UC, z*UC, Sens(im).db_image);
    hold on
    title(Sens(im).name)
    caxis([-db_range_for_output, 0])
    plot(probe_coords(:, 1)*UC, probe_coords(:, 3)*UC, 'go');
    for wall = 1:size(geometry, 1)
        plot(geometry(wall).coords(:, 1)*UC, geometry(wall).coords(:, 3)*UC, 'r')
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
filename_fig = sprintf('%s.fig', savename);
filename_mat = sprintf('%s.mat', savename);
saveas(fig, filename_fig)
% close all

time_4 = double(toc);

fprintf('Plotted in %.2f secs\n', time_4);
fprintf('%.2g MB used.\n', monitor_memory_whos);

times = [time_1, time_2, time_3, time_4];
% times.sens = Sens;
% times.views = Views;
% times.scat_info = image_block_info;

save(filename_mat, 'times', 'Sens', 'Views', 'image_block_info')

clear



end
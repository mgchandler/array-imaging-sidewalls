function fn_tfm(model_options)
% Computes the sensitivity maps for an array in contact with the solid
% block being inspected. Currently works for a rectangular block with some
% depth defined as the zmax location, and sidewall location defined as the
% xmax location.
%
%
% INPUTS:
% - model_options : struct (1, 1)
%       Structure containing options which the program will run with. For
%       more details on what each one is, see the fn_default_model_options
%       function. Possible options to include as fields are:
%       - pixel : double : DEFAULT = 5.0e-3
%       - probe_els : integer : DEFAULT = 32
%       - probe_angle : double : DEFAULT = 0
%       - probe_standoff : double : DEFAULT = 0
%       - probe_frequency : double : DEFAULT = 5.0e+6
%       - probe_pitch : double : DEFAULT = 1.0e-3
%       - el_length : double : DEFAULT = 0.0
%       - geom_shape : struct (1, 1)
%           - xmin : double : DEFAULT = -25.0e-3
%           - xmax : double : DEFAULT =  25.0e-3
%           - zmin : double : DEFAULT =   0.0e-3
%           - zmax : double : DEFAULT =  40.0e-3
%       - multi_freq : logical : DEFAULT = 0
%       - material_params : struct (1, 1)
%           - couplant_speed : double : DEFAULT = 340.0
%           - couplant_density : double : DEFAULT = 1.2
%           - solid_long_speed : double : DEFAULT = 6320.0
%           - solid_shear_speed : double : DEFAULT = 3130.0
%           - solid_density : double : DEFAULT = 2700.0
%       - scat_info : struct : DEFAULT sdh located at (0, 22e-3)
%       - boxsize : integer : DEFAULT = 0
%       - savepath : string : DEFAULT = ""
%       - savename : string : DEFAULT = "sens MODE GEOM VIEWS PITCH PIXEL WALLS"
%       - max_no_reflections : integer : DEFAULT = 1
%       - model_geometry : logical : DEFAULT = 0
%       - geometry : struct (no_walls, 1) : DEFAULT backwall



tic;

%% ---------------------------------------------------------------------- %
% Unpack model_config and model_options                                   %
% ---------------------------------------------------------------------- %%

PITCH = model_options.probe_pitch;
PIXEL = model_options.pixel;

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
scat_info = model_options.scat_info;
savepath = model_options.savepath;
savename = model_options.savename;
geometry = model_options.geometry;
max_num_reflections = model_options.max_no_reflections;
model_geometry = model_options.model_geometry;

% Additional parameters not directly dependent on inputs.
oversampling = 10;
no_cycles = 5;

no_walls = size(geometry, 1);

% Check whether we are in contact or immersion. If we are in contact, there
% will be no frontwall, and probe_standoff and probe_angle must equal zero.
% If we are in immersion, there must be a frontwall and probe_standoff must
% be non-zero.
is_frontwall = 0;
for wall = 1:no_walls
    if geometry(wall).name == "F"
        is_frontwall = 1;
        break
    end
end
if and(is_frontwall, probe_standoff ~= 0)
    is_contact = 0;
elseif and(~is_frontwall, and(probe_standoff == 0, probe_angle == 0))
    is_contact = 1;
else
    error('fn_sens: Invalid setup.')
end

xsize = xmax - xmin;
zsize = zmax - zmin;

UC = 1e3; % Unit conversion

clear PITCH is_frontwall wall



%% ------------------------------------------------------------------------
% Set up the scene.
% -------------------------------------------------------------------------

rot_matrix = [cos(probe_angle) 0 sin(probe_angle); 0 1 0; -sin(probe_angle) 0 cos(probe_angle)];

probe_coords = zeros(3, probe_els);
probe_coords(1, :) = linspace(0, (probe_els - 1) * el_length, probe_els);
probe_coords(1, :) = probe_coords(1, :) - mean(probe_coords(1, :));
probe_coords = probe_coords.' * rot_matrix;
probe_coords(:, 3) = probe_coords(:, 3) - probe_standoff;

clear rot_matrix


        
%% ------------------------------------------------------------------------
% Input signal
% -------------------------------------------------------------------------

% Create input signal.
time_step = 1 / (probe_frequency * oversampling); % What is the meaning of oversampling here?
max_t = 1.1 * 4 * sqrt(xsize ^ 2 + zsize ^ 2) / min(solid_long_speed, solid_shear_speed);
if probe_standoff ~= 0
    max_t = max_t + sqrt(probe_standoff^2 + (el_length*probe_els)^2) / couplant_speed;
end  
time_pts = ceil(max_t / time_step);
[~, ~, freq, in_freq_spec, fft_pts] = fn_create_input_signal(time_pts, probe_frequency, time_step , no_cycles);

clear probe_standoff oversampling no_cycles time_step max_t



%% ---------------------------------------------------------------------- %
% Scatterer Simulation path info - Contact                                %
% ---------------------------------------------------------------------- %%

tic;
    
mode_names = ["L", "T"];
speeds = [solid_long_speed, solid_shear_speed];

num_paths = 0;
for num_reflections_in_path = 0:max_num_reflections
    num_paths = num_paths + 2^(num_reflections_in_path + 1);
end

Path_info_list = repmat(fn_path_info( ...
    "", ...
    "", ...
    [0], ...
    0, ...
    0, ...
    [couplant_speed, solid_long_speed, solid_shear_speed], ...
    0, ...
    [1], ...
    [couplant_density, solid_density], ...
    0, ...
    0, ...
    probe_coords), ....
    num_paths, 1);

% If we are in the contact case
path = 1;
if is_contact
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
        Path_info_list(path) = Direct_path_info;
        path = path + 1;
        clear Direct_path_info
    end

    clear path_geometry
    
% ----------------------------------------------------------------------- %
% Skip Paths                                                              %
% ----------------------------------------------------------------------- %

    if max_num_reflections > 0
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
                    Path_info_list(path) = Skip_path_info;
                    path = path + 1;
                    clear Skip_path_info
                end
            end
            
            clear path_geometry
            
        end
    end

%% ---------------------------------------------------------------------- %
% Scatterer Simulation path info - Immersion                              %
% ---------------------------------------------------------------------- %%

% If we are in the immersion case.
elseif ~is_contact
    
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
        Path_info_list(path) = Direct_path_info;
        path = path + 1;
        clear Direct_path_info
    end
    
    clear path_geometry
    
% ----------------------------------------------------------------------- %
% Skip Paths                                                              %
% ----------------------------------------------------------------------- %

    if max_num_reflections > 0
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
                    Path_info_list(path) = Skip_path_info;
                    path = path + 1;
                    clear Skip_path_info
                end
            end
            
            clear path_geometry
            
        end
    end
    
end

time_1 = double(toc);

fprintf('Setup time %.2f secs.\n', time_1);

clear max_num_reflections no_walls mode_names speeds num_reflections_in_path
clear path mode_name wall mode1 mode1_name mode2 mode2_name


%% ---------------------------------------------------------------------- %
% Create simulation views                                                 %
% ---------------------------------------------------------------------- %%
% Precompute views for all scatterers, as well as imaging views.

tic;

% Compute Imaging Paths
Paths = repmat(fn_compute_ray(scat_info, Path_info_list(1), probe_frequency), 1, num_paths);
for path = 2:num_paths
    Paths(path) = fn_compute_ray(scat_info, Path_info_list(path), probe_frequency);
end

% Create views from these paths.
Views = fn_make_views(Paths, 0);

clear Paths



%% ---------------------------------------------------------------------- %
% Geometry setup.   NEEDS WORK                                            %
% ---------------------------------------------------------------------- %%

% If we are in the contact case. There is no frontwall.
if is_contact
    if model_geometry
        backwall_views = fn_backwall_views( ...
            probe_coords, probe_angle, 0, back_wall, [couplant_speed solid_long_speed solid_shear_speed], ...
            [couplant_density solid_density], probe_frequency, el_length ...
        );
    end

% If we are in the immersion case. There will be a frontwall reflection.
else
    if model_geometry
        frontwall_view = fn_frontwall_view( ...
            probe_coords, probe_angle, front_wall, [couplant_speed solid_long_speed solid_shear_speed], ...
            [couplant_density solid_density], probe_frequency, el_length ...
        );
    
        backwall_views = fn_backwall_views( ...
            probe_coords, probe_angle, front_wall, back_wall, [couplant_speed solid_long_speed solid_shear_speed], ...
            [couplant_density solid_density], probe_frequency, el_length ...
        );
    end
    
end

time_2 = double(toc);

fprintf('Rays traced in %.2f secs.\n', time_2);

clear probe_angle probe_frequency el_length couplant_speed couplant_density
clear solid_long_speed solid_shear_speed solid_density



%% ---------------------------------------------------------------------- %
% Simulation.                                                             %
% ---------------------------------------------------------------------- %%

tic

% Scatterer simulation.
out_freq_spec = 0;
num_scatterers = size(scat_info.image_block, 1);

for scatterer = 1 : num_scatterers
    for view = 1 : size(Views, 1)
        weights = Views(view).weights(:, scatterer, 1);
        scat_amp = Views(view).scat_amps(:, scatterer, 1);
        amp = conj(scat_amp .* weights);
        out_freq_spec = out_freq_spec + ...
            fn_propagate_spectrum_mc(freq, in_freq_spec, Views(view).min_times(:, scatterer), amp, 0);
        
        clear weights scat_amp amp
    end
end

if model_geometry
    for view = 1 : size(backwall_views, 1)
        out_freq_spec = out_freq_spec + ...
            fn_propagate_spectrum_mc(freq, in_freq_spec, backwall_views(view).min_times, backwall_views(view).weights, 0);
    end
    % If immersion, then do frontwall.
    if ~is_contact
        out_freq_spec = out_freq_spec + ...
            fn_propagate_spectrum_mc(freq, in_freq_spec, frontwall_view.min_times, frontwall_view.weights, 0);
        
        clear frontwall_view
    end
    clear backwall_views
end

% Convert back to time.
[FMC_time, FMC_time_data] = fn_convert_spectrum_to_time(freq, out_freq_spec, fft_pts, time_pts);
FMC_time = FMC_time';

% Hilbert Filtering (?) from fast_DAS function
diagonals = spdiags([1:length(FMC_time)]' < length(FMC_time)/2, 0, length(FMC_time), length(FMC_time));
FMC_time_data = ifft(diagonals * fft(FMC_time_data));

time_3 = double(toc);

fprintf('Simulated in in %.2f secs\n', time_3);

clear model_geometry time_pts freq in_freq_spec fft_pts is_contact out_freq_spec
clear num_scatterers scatterer view diagonals



%% ---------------------------------------------------------------------- %
% Imaging Setup                                                           %
% ---------------------------------------------------------------------- %%

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
scatterer_coords = reshape(image_block, zpts+1, xpts+1, 3);
are_points_in_geometry = (scatterer_coords(:,:,1) >= xmin) .* ...
                         (scatterer_coords(:,:,1) <= xmax) .* ...
                         (scatterer_coords(:,:,3) >= zmin) .* ...
                         (scatterer_coords(:,:,3) < zmax);

Paths_im = repmat(fn_compute_ray(image_block_info, Path_info_list(1)), 1, num_paths);
for path = 2:num_paths
    Paths_im(path) = fn_compute_ray(image_block_info, Path_info_list(path));
end

Views_im = fn_make_views(Paths_im, 1);

Number_of_ims = size(Views_im, 1);
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
    error('fn_sens: Unexpected number of images being plotted.\n%d image(s) being plotted.', Number_of_ims)
end

Ims = repmat(fn_create_im("-", xpts+1, zpts+1), Number_of_ims, 1);
for view = 1 : Number_of_ims
    Ims(view).name = Views(view).name;
end

clear PIXEL xmin xmax zmin zmax xsize zsize num_paths Path_info_list image_block
clear image_block_info scatterer_coords Paths_im



%% ---------------------------------------------------------------------- %
% Imaging                                                                 %
% ---------------------------------------------------------------------- %%

% Lookup times in FMC data.
for tr_pair = 1 : probe_els ^ 2
    for view = 1 : Number_of_ims
        tau = reshape(Views_im(view).min_times(tr_pair, :), [zpts+1, xpts+1]);
        Ims(view).image = Ims(view).image + are_points_in_geometry .* ...
            interp1(FMC_time, FMC_time_data(:, tr_pair), tau, 'linear', 0);
    end
end

time_4 = double(toc);

fprintf('Imaged in %.2f secs\n', time_4);

clear FMC_time FMC_time_data xpts zpts are_points_in_geometry tr_pair view tau



%% ---------------------------------------------------------------------- %
% Plotting.                                                               %
% ---------------------------------------------------------------------- %%

tic

max_ = 0;
for view = 1 : Number_of_ims
    if max(abs(Ims(view).image(:))) > max_
        max_ = max(abs(Ims(view).image(:)));
    end
end

for view = 1 : Number_of_ims
   Ims(view).db_image = 20 * log10(abs(Ims(view).image) ./ max_); 
end

% Plot.
fig = figure(1);
ax = repmat(subplot(plot_z, plot_x, 1), Number_of_ims, 1);
for im = 1:Number_of_ims
    ax(im) = subplot(plot_z, plot_x, im);
    imagesc(im_x*UC, im_z*UC, Ims(im).db_image);
    hold on
    title(Ims(im).name)
    caxis([-db_range_for_output, 0])
    plot(probe_coords(:, 1)*UC, probe_coords(:, 3)*UC, 'go');
    for wall = 1:size(geometry, 1)
        plot(geometry(wall).coords(:, 1)*UC, geometry(wall).coords(:, 3)*UC, 'r')
    end
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



time_5 = double(toc);

fprintf('Plotted in %.2f secs\n', time_5);

times = [time_1, time_2, time_3, time_4, time_5];



if savepath ~= ""
    cd(savepath)
    filename_fig = sprintf('%s.fig', savename);
    filename_mat = sprintf('%s.mat', savename);
    savefig(filename_fig)
    save(filename_mat, 'times', 'Ims', 'Views')
end
% close all

clear



end
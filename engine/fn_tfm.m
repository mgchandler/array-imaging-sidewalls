function fn_tfm(model_options)
% Computes the sensitivity maps for an array in contact with the solid
% block being inspected. Currently works for a rectangular block with some
% depth defined as the zmax location, and sidewall location defined as the
% xmax location.
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
xmin = model_options.geom_shape.xmin;
xmax = model_options.geom_shape.xmax;
zmin = model_options.geom_shape.zmin;
zmax = model_options.geom_shape.zmax;
scat_info = model_options.scat_info;
savepath = model_options.savepath;
savename = model_options.savename;
geometry = model_options.geometry;
wall_for_imaging = model_options.wall_for_imaging;

probe_angle = model_options.probe_angle;
probe_standoff = model_options.probe_standoff;
probe_frequency = model_options.probe_frequency;
el_length = model_options.el_length;
couplant_speed = model_options.material_params.couplant_speed;
couplant_density = model_options.material_params.couplant_density;
solid_long_speed = model_options.material_params.solid_long_speed;
solid_shear_speed = model_options.material_params.solid_shear_speed;
solid_density = model_options.material_params.solid_density;
max_num_reflections = model_options.max_no_reflections;
model_geometry = model_options.model_geometry;
multi_freq = model_options.multi_freq;
norm_to = model_options.norm_to;
db_range_for_output = model_options.db_range_for_output;

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

clear is_frontwall wall



%% ------------------------------------------------------------------------
% Set up the scene.
% -------------------------------------------------------------------------

rot_matrix = [cos(probe_angle) 0 sin(probe_angle); 0 1 0; -sin(probe_angle) 0 cos(probe_angle)];

probe_coords = zeros(3, probe_els);
probe_coords(1, :) = linspace(0, (probe_els - 1) * PITCH, probe_els);
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
    max_t = max_t + sqrt(probe_standoff^2 + (PITCH*probe_els)^2) / couplant_speed;
end  
time_pts = ceil(max_t / time_step)+1;
[~, ~, freq, in_freq_spec, fft_pts] = fn_create_input_signal(time_pts, probe_frequency, time_step , no_cycles);

clear probe_standoff oversampling no_cycles time_step



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
    [1], ...
    [couplant_density, solid_density], ...
    0, ...
    0, ...
    0, ...
    probe_coords), ...
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
            [1], ... % index for material identity
            [couplant_density, solid_density], ...
            probe_frequency, ...
            PITCH, ...
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
...%                %%  Use these names to include the wall in the path name
                        sprintf("%s %s %s", mode1_name, path_geometry.name, mode2_name), ...
                        sprintf("%s %s %s", mode2_name, path_geometry.name, mode1_name), ...
...%                %%  Use these names to exclude the wall in the path name
...%                         sprintf("%s %s", mode1_name, mode2_name), ...
...%                         sprintf("%s %s", mode2_name, mode1_name), ...
                        [mode1, mode2], ...
                        path_geometry, ...
                        [speeds(mode1+1), speeds(mode2+1)], ...
                        [couplant_speed, solid_long_speed, solid_shear_speed], ...
                        [1, 1], ...
                        [couplant_density, solid_density], ...
                        probe_frequency, ...
                        PITCH, ...
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

    wall_names = repmat("", size(geometry, 1), 1);
    for wall = 1:size(geometry, 1)
        wall_names(wall) = geometry(wall).name;
    end
    where_F = logical(count(wall_names, "F"));
    path_geometry = geometry(where_F);

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
            [0, 1], ...
            [couplant_density, solid_density], ...
            probe_frequency, ...
            PITCH, ...
            el_length, ...
            probe_coords ...
        );
        Path_info_list(path) = Direct_path_info;
        path = path + 1;
        clear Direct_path_info
    end

    clear path_geometry wall_names

% ----------------------------------------------------------------------- %
% Skip Paths                                                              %
% ----------------------------------------------------------------------- %

    if max_num_reflections > 0
        non_fw_idxs = [1:no_walls];
        non_fw_idxs = non_fw_idxs(~where_F);
        for wall2 = 1:no_walls-1
            wall = non_fw_idxs(wall2);

            path_geometry = repmat(geometry(where_F), 2, 1);
            path_geometry(2) = geometry(wall);
            for mode1 = 0:1 % Mode of the first leg
                mode1_name = mode_names(mode1+1);
                for mode2 = 0:1 % Mode of the second leg
                    mode2_name = mode_names(mode2+1);
                    Skip_path_info = fn_path_info( ...
...%                %%  Use these names to include the wall in the path name
                        sprintf("%s %s %s", mode1_name, path_geometry(2).name, mode2_name), ...
                        sprintf("%s %s %s", mode2_name, path_geometry(2).name, mode1_name), ...
...%                %%  Use these names to exclude the wall in the path name
...%                         sprintf("%s %s", mode1_name, mode2_name), ...
...%                         sprintf("%s %s", mode2_name, mode1_name), ...
                        [0, mode1, mode2], ...
                        path_geometry, ...
                        [couplant_speed, speeds(mode1+1), speeds(mode2+1)], ...
                        [couplant_speed, solid_long_speed, solid_shear_speed], ...
                        [0, 1, 1], ...
                        [couplant_density, solid_density], ...
                        probe_frequency, ...
                        PITCH, ...
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

fn_print_time('Setup', time_1)

clear max_num_reflections no_walls mode_names speeds num_reflections_in_path
clear path mode_name wall mode1 mode1_name mode2 mode2_name wall2
    

% If no FMC data supplied, then skip simulation.
if ~isstruct(model_options.FMC_data)


    %% ---------------------------------------------------------------------- %
    % Create simulation views                                                 %
    % ---------------------------------------------------------------------- %%
    % Precompute views for all scatterers, as well as imaging views.

    tic;

    % Compute Imaging Paths
    if multi_freq
        frequency = freq(2:end);
    else
        frequency = probe_frequency;
    end
    Paths = repmat(fn_compute_ray(scat_info, Path_info_list(1), geometry, frequency), 1, num_paths);
    path = 1;
    ii = 1;
    while path < num_paths
        ii = ii + 1;
        if length(Path_info_list(ii).speeds) > 1
            if wall_for_imaging == Path_info_list(ii).path_geometry.name
                path = path + 1;
                Paths(path) = fn_compute_ray(scat_info, Path_info_list(ii), geometry, frequency);
            end
        else % Must be direct
            path = path+1;
            Paths(path) = fn_compute_ray(scat_info, Path_info_list(ii), geometry, frequency);
        end
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

            backwall_views = fn_make_geometry_views(probe_coords, geometry, ...
                [couplant_speed solid_long_speed solid_shear_speed], ...
                [couplant_density solid_density], probe_frequency, PITCH, el_length, 1 ...
            );
        end

    % If we are in the immersion case. There will be a frontwall reflection.
    else
        if model_geometry

            backwall_views = fn_make_geometry_views(probe_coords, geometry, ...
                [couplant_speed solid_long_speed solid_shear_speed], ...
                [couplant_density solid_density], probe_frequency, PITCH, el_length, 1 ...
            );
        end
    end

    time_2 = double(toc);
    
    fn_print_time('Rays traced', time_2)

    clear probe_angle probe_frequency el_length couplant_speed couplant_density
    clear solid_long_speed solid_shear_speed solid_density PITCH



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
            valid_path = Views(view).valid_path(:, scatterer);
            amp = (scat_amp .* weights .* valid_path);
            out_freq_spec = out_freq_spec + ...
                fn_propagate_spectrum_mc(freq, in_freq_spec, Views(view).min_times(:, scatterer), amp, 0);

            clear weights scat_amp amp
        end
    end

    if model_geometry
        [num_geom_views, ~] = size(backwall_views);
        for view = 1 : num_geom_views
            bw_amp = conj(backwall_views(view).weights .* backwall_views(view).valid_path);
            out_freq_spec = out_freq_spec + ...
                fn_propagate_spectrum_mc(freq, in_freq_spec, backwall_views(view).min_times, bw_amp, 0);
        end
        clear backwall_views
    end

    % Convert back to time.
    [FMC_time, FMC_time_data] = fn_convert_spectrum_to_time(freq, out_freq_spec, fft_pts, time_pts);
    FMC_time = FMC_time';

    % Hilbert Filtering (?) from fast_DAS function
    diagonals = spdiags([1:length(FMC_time)]' < length(FMC_time)/2, 0, length(FMC_time), length(FMC_time));
    FMC_time_data = ifft(diagonals * fft(FMC_time_data));
    
    plot_time = linspace(0, max_t, time_pts)*10^6;
    figure(3)
    imagesc(plot_time, [1:1024], abs(FMC_time_data)')
    xlabel('Time (us)')
    ylabel('tx-rx Index')
    savefig(sprintf('%s_FMC.fig', savename));
    xsize = xmax - xmin;
    zsize = zmax - zmin;

    time_3 = double(toc);
    
    fn_print_time('Simulated', time_3)

    clear model_geometry time_pts freq in_freq_spec fft_pts out_freq_spec
    clear num_scatterers scatterer view diagonals
    
    
    
%% ---------------------------------------------------------------------- %
% Load FMC data                                                           %
% ---------------------------------------------------------------------- %%
else
    xsize = xmax - xmin;
    zsize = zmax - zmin;
    
    FMC_time = model_options.FMC_data.time;
    FMC_time_data = model_options.FMC_data.data;
    
    figure(2)
    FMC_for_plotting = abs(FMC_time_data');
    FMC_for_plotting(:, 1:751) = 0;
    imagesc(FMC_time(251:end)*10^6, [1:32], FMC_for_plotting)
    xlabel('Time (us)')
    ylabel('tx-rx Index')
end



%% ---------------------------------------------------------------------- %
% Imaging Setup                                                           %
% ---------------------------------------------------------------------- %%

tic

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
thispath = 1;
for path = 2:size(Path_info_list, 1)
    % Only get the paths we want to image. This will always include direct
    % paths (when path_geometry field == 0), and may include some extra
    % ones if we are modelling wall reflections.
    if is_contact
        if ~isstruct(Path_info_list(path).path_geometry) % If direct contact
            Paths_im(thispath) = fn_compute_ray(image_block_info, Path_info_list(path));
            thispath = thispath + 1;
        elseif Path_info_list(path).path_geometry(end).name == wall_for_imaging
            Paths_im(thispath) = fn_compute_ray(image_block_info, Path_info_list(path));
            thispath = thispath + 1;
        end
    else
        if Path_info_list(path).path_geometry(end).name == "F" % If direct immersion
            Paths_im(thispath) = fn_compute_ray(image_block_info, Path_info_list(path));
            thispath = thispath + 1;
        elseif Path_info_list(path).path_geometry(end).name == wall_for_imaging
            Paths_im(thispath) = fn_compute_ray(image_block_info, Path_info_list(path));
            thispath = thispath + 1;
        end
    end
end

Views_im = fn_make_views(Paths_im, 1);

Number_of_ims = size(Views_im, 1);
if Number_of_ims == 3
    plot_x = 3;
    plot_z = 1;
elseif Number_of_ims == 21
    plot_x = 3;
    plot_z = 7;
elseif Number_of_ims == 55
    plot_x = 11;
    plot_z = 5;
else
    error('fn_sens: Unexpected number of images being plotted.\n%d image(s) being plotted.', Number_of_ims)
end

Ims = repmat(fn_create_im("-", xpts+1, zpts+1), Number_of_ims, 1);
for view = 1 : Number_of_ims
    Ims(view).name = Views_im(view).name;
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

fn_print_time('Imaged', time_4)

clear FMC_time FMC_time_data xpts zpts are_points_in_geometry tr_pair view tau



%% ---------------------------------------------------------------------- %
% Plotting.                                                               %
% ---------------------------------------------------------------------- %%

tic

if norm_to == 0
    max_ = 0;
    for view = 1 : Number_of_ims
        if max(abs(Ims(view).image(:))) > max_
            max_ = max(abs(Ims(view).image(:)));
        end
    end
else
    max_ = norm_to;
end

for view = 1 : Number_of_ims
   Ims(view).db_image = 20 * log10(abs(Ims(view).image) ./ max_); 
end

% Plot.
fig = figure(1);
% ax = repmat(subplot(plot_z, plot_x, 1), Number_of_ims, 1);
% sgtitle(sprintf('p = %.2f', (probe_coords(2,1)-probe_coords(1,1))*10^3))
t = tiledlayout(plot_z, plot_x, 'TileSpacing', 'Compact');
for im = 1:Number_of_ims
    h(im) = nexttile;
%     ax(im) = subplot(plot_z, plot_x, im);
    imagesc(im_x*UC, im_z*UC, Ims(im).db_image);
    hold on
    title(Ims(im).name)
    caxis([-db_range_for_output, 0])
    plot(probe_coords(:, 1)*UC, probe_coords(:, 3)*UC, 'go');
    for wall = 1:size(geometry, 1)
        plot(geometry(wall).coords(:, 1)*UC, geometry(wall).coords(:, 3)*UC, 'r')
    end
%     for s = 1 : size(scat_info.image_block, 1)
%         if scat_info.type ~= 'image'
%             rectangle('Position', [scat_info.image_block(s, 1)*UC - 1, scat_info.image_block(s, 3)*UC - 1, 2, 2], 'EdgeColor', 'r');
%         end
%     end
    
    if mod(im, plot_x) ~= 1
        set(gca, 'yticklabel', {[]})
    end
    if im <= Number_of_ims - plot_x
        set(gca, 'xticklabel', {[]})
    end
    
    axis equal; axis tight;
end
xlabel(t, 'x (mm)')
ylabel(t, 'z (mm)')

c = colorbar(h(1), 'AxisLocation','in');
c.Layout.Tile = 'north';
c.Label.String = 'dB';

% % Resize image and move for colorbar
% set(fig, 'Position', [20, 20, im_width, im_height])
% posa = cell2mat(get(ax, 'Position'));
% h = colorbar;
% for im = 1 : Number_of_ims
%     set(ax(im), 'Position', [posa(im, 1), posa(im, 2)*0.8, posa(im, 3)*1.1, posa(im, 4)*1.1])
% end
% set(ax, 'units', 'pix')
% set(h, 'units', 'pix')
% posf = get(fig, 'Position'); % gives x left, y bottom, width, height
% set(fig, 'Position',  [posf(1:2) posf(3)*1.1 posf(4)])
% hpos = h.Position;
% posa = cell2mat(get(ax, 'Position'));
% pos1 = posa(1, :);
% set(h, 'Position', [hpos(1)+10, hpos(2), hpos(3)*2, pos1(2)-hpos(2)+pos1(4)])
% 
% h.Label.String = 'dB';



time_5 = double(toc);

fn_print_time('Plotted', time_5)

if isstruct(model_options.FMC_data)
    times = [time_1];
else
    times = [time_1, time_2, time_3, time_4, time_5];
end



if savepath ~= ""
    cd(savepath)
    filename_fig = sprintf('%s.fig', savename);
    filename_mat = sprintf('%s.mat', savename);
    savefig(filename_fig)
%     if exist('Views')
%         save(filename_mat, 'times', 'Ims', 'Views')
%     else
    save(filename_mat, 'times', 'Ims')
%     end
end
% close all

clear



end
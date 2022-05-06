function fn_tfm(model_options)
% Computes the sensitivity maps for an array in contact with the solid
% block being inspected. Currently works for a rectangular block with some
% depth defined as the zmax location, and sidewall location defined as the
% xmax location.
%
% INPUTS:
% - model_options : struct 
%       Structure containing options which will overwrite the defaults.
%       Possible options to include as fields are:
%       - data : struct : DEFAULT none
%           If FMC data generated externally, it can be passed in using
%           this option. The simulation will be skipped and the program
%           will proceed straight to the TFM.
%       - material : struct
%           - couplant_density : double : DEFAULT 1.2
%           - couplant_v : double : DEFAULT 340.0
%           - density : double : DEFAULT 2700.0
%           - modulus : double : DEFAULT 70e9
%           - poisson : double : DEFAULT 0.34
%       - mesh : struct
%           - geom : struct
%               - x : double array
%                   List of x coordinates
%               - y : double array
%                   List of z coordinates
%               - geometry : struct
%                   Output from fn_make_geometry()
%           - scat : struct
%               Output from fn_scat_info()
%       - model : struct
%           - boxsize : double : DEFAULT 0.0
%           - db_range : double : DEFAULT 40.0
%           - max_no_reflections : integer : DEFAULT 1
%           - model_geom : logical : DEFAULT 1
%           - multi_freq : logical : DEFAULT 0
%           - norm_to : double : DEFAULT 0
%           - npw : double : DEFAULT 45
%           - pixel : double : DEFAULT 0.5e-3
%           - savename : string : DEFAULT "TFM-Sens Image Plot"
%           - savepath : string : DEFAULT ""
%           - wall_for_imaging : string : DEFAULT "B1"
%       - probe : struct
%           - angle : double : DEFAULT 0
%           - freq : double : DEFAULT 5.0e+6
%           - num_els : integer : DEFAULT 32
%           - standoff : double : DEFAULT 0
%           - separation : double : DEFAULT 0.05e-3
%           - width : double : DEFAULT 0.45e-3



tic;

%% ---------------------------------------------------------------------- %
% Unpack model_config and model_options                                   %
% ---------------------------------------------------------------------- %%

PITCH = model_options.probe.width + model_options.probe.separation;
PIXEL = model_options.model.pixel;

probe_els = model_options.probe.num_els;
xmin = min(cell2mat(model_options.mesh.geom.x));
xmax = max(cell2mat(model_options.mesh.geom.x));
zmin = min(cell2mat(model_options.mesh.geom.z));
zmax = max(cell2mat(model_options.mesh.geom.z));
scat_info = model_options.mesh.scat;
savepath = model_options.model.savepath;
savename = model_options.model.savename;
geometry = model_options.mesh.geom.geometry;
wall_for_imaging = model_options.model.wall_for_imaging;
boxsize = model_options.model.boxsize;

probe_angle = model_options.probe.angle;
probe_standoff = model_options.probe.standoff;
probe_frequency = model_options.probe.freq;
el_length = model_options.probe.width;
couplant_speed = model_options.material.couplant_v;
couplant_density = model_options.material.couplant_density;
solid_long_speed = sqrt(model_options.material.modulus * (1 - model_options.material.poisson) / (model_options.material.density * (1 + model_options.material.poisson) * (1 - 2*model_options.material.poisson)));
solid_shear_speed = sqrt(model_options.material.modulus / (2 * model_options.material.density * (1 + model_options.material.poisson)));
solid_density = model_options.material.density;
max_num_reflections = model_options.model.max_no_reflections;
norm_to = model_options.model.norm_to;
db_range_for_output = model_options.model.db_range;
npw = model_options.mesh.n_per_wl;

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



%% ------------------------------------------------------------------------
% Get probe for plotting later
% -------------------------------------------------------------------------

rot_matrix = [cos(probe_angle) 0 sin(probe_angle); 0 1 0; -sin(probe_angle) 0 cos(probe_angle)];

probe_coords = zeros(3, probe_els);
probe_coords(1, :) = linspace(0, (probe_els - 1) * PITCH, probe_els);
probe_coords(1, :) = probe_coords(1, :) - mean(probe_coords(1, :));
probe_coords = probe_coords.' * rot_matrix;
probe_coords(:, 3) = probe_coords(:, 3) - probe_standoff;

clear rot_matrix



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
    probe_coords, ...
    npw), ...
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
            probe_coords, ...
            npw ...
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
                        probe_coords, ...
                        npw ...
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
            probe_coords, ...
            npw ...
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
                        probe_coords, ...
                        npw ...
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
    


%% ---------------------------------------------------------------------- %
% Simulate FMC data                                                       %
% ---------------------------------------------------------------------- %%
% If no FMC data supplied, then skip simulation.
if ~isstruct(model_options.data)

    [FMC_time, FMC_time_data] = fn_simulate_fmc(model_options);

%% ---------------------------------------------------------------------- %
% Load FMC data                                                           %
% ---------------------------------------------------------------------- %%
else
    xsize = xmax - xmin;
    zsize = zmax - zmin;
    
    FMC_time = model_options.data.time;
    FMC_time_data = model_options.data.data;
    
    FMC_for_plotting = abs(FMC_time_data);
    FMC_for_plotting(1:751, :) = 0;
%     fn_plot_FMC_at_time(FMC_for_plotting, FMC_time, Path_info_list(3), Path_info_list(6), [scat_info.x, scat_info.y, scat_info.z], sprintf('%s_FMC.png', savename));
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

image_block_info = fn_scat_info("image", image_block(:, 1), image_block(:, 2), image_block(:, 3));
scatterer_coords = reshape(image_block, zpts+1, xpts+1, 3);
are_points_in_geometry = (scatterer_coords(:,:,1) >= xmin) .* ...
                         (scatterer_coords(:,:,1) <= xmax) .* ...
                         (scatterer_coords(:,:,3) >= zmin) .* ...
                         (scatterer_coords(:,:,3) <  zmax);

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

clear PIXEL xsize zsize num_paths Path_info_list image_block
clear image_block_info scatterer_coords Paths_im



%% ---------------------------------------------------------------------- %
% Imaging                                                                 %
% ---------------------------------------------------------------------- %%

% Lookup times in FMC data.

if strcmp(model_options.model.interp_method, 'lanczos')
    for tr_pair = 1 : probe_els ^ 2
        for view = 1 : Number_of_ims
            tau = reshape(Views_im(view).min_times(tr_pair, :), [zpts+1, xpts+1]);
            Ims(view).image = Ims(view).image + are_points_in_geometry .* ...
                fn_lanczos_interp(FMC_time, FMC_time_data(:, tr_pair), tau, 3);
        end
    end
else % Assume linear interpolation as default.
    for tr_pair = 1 : probe_els ^ 2
        for view = 1 : Number_of_ims
            tau = reshape(Views_im(view).min_times(tr_pair, :), [zpts+1, xpts+1]);
            Ims(view).image = Ims(view).image + are_points_in_geometry .* ...
                interp1(FMC_time, FMC_time_data(:, tr_pair), tau, 'linear', 0);
        end
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

% im_idxs = [4, 12, 17];

% Plot.
fig = figure(1);
t = tiledlayout(plot_z, plot_x, 'TileSpacing', 'Compact');
for im = 1:Number_of_ims
%     im = im_idxs(im1);
    h(im) = nexttile;
    imagesc(im_x*UC, im_z*UC, Ims(im).db_image);
    hold on
    title(Ims(im).name)
    caxis([-db_range_for_output, 0])
    plot(probe_coords(:, 1)*UC, probe_coords(:, 3)*UC, 'go');
    for wall = 1:size(geometry, 1)
        plot(geometry(wall).coords(:, 1)*UC, geometry(wall).coords(:, 3)*UC, 'r')
    end
    if boxsize ~= 0
        for s = 1 : size(scat_info.x, 1)
            if strcmp(scat_info.type, 'image')
                rectangle('Position', [scat_info.x(s)*UC - boxsize*UC/2, scat_info.z(s)*UC - boxsize*UC/2, boxsize*UC, boxsize*UC], 'EdgeColor', 'r');
            end
        end
    end
    
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



if savepath ~= ""
    cd(savepath)
    filename_fig = sprintf('%s.fig', savename);
    filename_mat = sprintf('%s.mat', savename);
    savefig(filename_fig)
    save(filename_mat, 'Ims')
end

clear



end
function fn_sens(model_options)
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
%               Output from fn_scat_info() function.
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
% Mins and maxes used for grid, so make it slightly smaller than the geometry.
xmin = min(cell2mat(model_options.mesh.geom.x)) + 0.1e-3;
xmax = max(cell2mat(model_options.mesh.geom.x)) - 0.1e-5;
zmin = min(cell2mat(model_options.mesh.geom.z)) + 0.1e-5;
zmax = max(cell2mat(model_options.mesh.geom.z)) - 0.1e-5;
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
solid_long_speed = 6320;%sqrt(model_options.material.modulus * (1 - model_options.material.poisson) / (model_options.material.density * (1 + model_options.material.poisson) * (1 - 2*model_options.material.poisson)));
solid_shear_speed = 3130;%sqrt(model_options.material.modulus / (2 * model_options.material.density * (1 + model_options.material.poisson)));
solid_density = model_options.material.density;
max_num_reflections = model_options.model.max_no_reflections;
model_geometry = model_options.model.model_geom;
multi_freq = model_options.model.multi_freq;
norm_to = model_options.model.norm_to;
db_range_for_output = model_options.model.db_range;
npw = model_options.mesh.n_per_wl;
image_block_info = model_options.mesh.scat;

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

if boxsize ~= 0
    warning('fn_sens: boxsize=%0.2g is ignored', boxsize)
end

xsize = xmax - xmin;
zsize = zmax - zmin;

UC = 1e3; % Unit conversion

clear is_frontwall wall



%% ---------------------------------------------------------------------- %
% Set up the scene.                                                       %
% ---------------------------------------------------------------------- %%

rot_matrix = [cos(probe_angle) 0 sin(probe_angle); 0 1 0; -sin(probe_angle) 0 cos(probe_angle)];

probe_coords = zeros(3, probe_els);
probe_coords(1, :) = linspace(0, (probe_els - 1) * el_length, probe_els);
probe_coords(1, :) = probe_coords(1, :) - mean(probe_coords(1, :));
probe_coords = probe_coords.' * rot_matrix;
probe_coords(:, 3) = probe_coords(:, 3) - probe_standoff - 1e-5;

clear probe_angle rot_matrix



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

clear num_reflections_in_path

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
    0), ....
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
            % For brief testing, only want reflections from the big
            % sidewall.
            if geometry(wall).name ~= wall_for_imaging
                continue
            end
            path_geometry = geometry(wall);
            for mode1 = 0:1 % Mode of the first leg
                mode1_name = mode_names(mode1+1);
                for mode2 = 0:1 % Mode of the second leg
                    mode2_name = mode_names(mode2+1);
                    Skip_path_info = fn_path_info( ...
...%                         sprintf("%s %s %s", mode1_name, path_geometry.name, mode2_name), ...
...%                         sprintf("%s %s %s", mode2_name, path_geometry.name, mode1_name), ...
                        sprintf("%s %s", mode1_name, mode2_name), ...
                        sprintf("%s %s", mode2_name, mode1_name), ...
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
    
% ----------------------------------------------------------------------- %
% Skip Paths                                                              %
% ----------------------------------------------------------------------- %

    if max_num_reflections > 0
        for wall = 2:no_walls
            if geometry(wall).name ~= wall_for_imaging
                continue
            end
            path_geometry = repmat(geometry(1), 2, 1);
            path_geometry(2) = geometry(wall);
            for mode1 = 0:1 % Mode of the first leg
                mode1_name = mode_names(mode1+1);
                for mode2 = 0:1 % Mode of the second leg
                    mode2_name = mode_names(mode2+1);
                    Skip_path_info = fn_path_info( ...
...%                         sprintf("%s %s %s", mode1_name, path_geometry(2).name, mode2_name), ...
...%                         sprintf("%s %s %s", mode2_name, path_geometry(2).name, mode1_name), ...
                        sprintf("%s %s", mode1_name, mode2_name), ...
                        sprintf("%s %s", mode2_name, mode1_name), ...
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
        end
    end
    
end

time_1 = double(toc);

fn_print_time('Setup', time_1)

clear el_length couplant_speed couplant_density solid_long_speed solid_shear_speed solid_density
clear max_num_reflections no_walls is_contact mode_names speeds mode mode_name wall mode1 mode2
clear mode1_name mode2_name PITCH



%% ---------------------------------------------------------------------- %
% Scatterer and Imaging setup                                             %
% ---------------------------------------------------------------------- %%
% Precompute views for all scatterers, as well as imaging views.

tic;

db_range_for_output = 40;

% The total number of points where scatterers will be placed in the x- and
% z- axes.
xpts = round(xsize / PIXEL);
zpts = round(zsize / PIXEL);
% Helper arrays for use in the meshgrid.
x = linspace(xmin, xmax, xpts+1);
z = linspace(zmin, zmax, zpts+1);
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
image_block_info.x = image_block(:, 1);
image_block_info.y = image_block(:, 2);
image_block_info.z = image_block(:, 3);

% Compute Imaging Paths
Paths = repmat(fn_compute_ray(image_block_info, Path_info_list(1), geometry, probe_frequency), 1, num_paths);
for path = 2:num_paths
    Paths(path) = fn_compute_ray(image_block_info, Path_info_list(path), geometry, probe_frequency);
end
clear probe_frequency num_paths path Path_info_list boxpix X Z pt xpt zpt path

% Create views from these paths.
Views = fn_make_views(Paths, 1);
Number_of_ims = size(Views, 1);
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

Sens = repmat(fn_create_im("-", xpts+1, zpts+1), Number_of_ims, 1);
for view = 1 : Number_of_ims
    Sens(view).name = Views(view).name;
end

time_2 = double(toc);

fn_print_time('Rays traced', time_2)

clear Paths



%% ---------------------------------------------------------------------- %
% Sensitivity Images                                                      %
% ---------------------------------------------------------------------- %%

tic

% Obscure points outside of the geometry: these will be zero in logical
% checks below, so multiply out all checks to kill pixels introduced by
% boxsize variable.
% are_points_in_geometry = (scatterer_coords(:,:,1) >= xmin) .* ...
%                          (scatterer_coords(:,:,1) <= xmax) .* ...
%                          (scatterer_coords(:,:,3) >= zmin) .* ...
%                          (scatterer_coords(:,:,3) < zmax);
are_points_in_geometry = (image_block(:,1) >= xmin) .* ...
                         (image_block(:,1) <= xmax) .* ...
                         (image_block(:,3) >= zmin) .* ...
                         (image_block(:,3) < zmax);

clear xmin xmax zmin zmax scatterer_coords



%% ---------------------------------------------------------------------- %
% Start of sensitivity loop                                               %
% ---------------------------------------------------------------------- %%

tic;

% ----------------------------------------------------------------------- %
% Simulation Step                                                         %
% ----------------------------------------------------------------------- %

for view = 1 : Number_of_ims
    weights = Views(view).weights;
    scat_amp = Views(view).scat_amps;
    valid_path = Views(view).valid_path;
    Sens(view).image = reshape(sum(conj(scat_amp .* weights .* valid_path .* are_points_in_geometry'), 1)/size(weights, 1), zpts+1, xpts+1);
end

time_3 = double(toc);

fn_print_time('Finished Sensitivity Loop', time_3)

clear PIXEL probe_els boxsize xsize zsize time_pts freq in_freq_spec fft_pts
clear xpts zpts are_points_in_geometry sens_i_min sens_i_max sens_k_min sens_k_max
clear grid_pt xpt_im zpt_im view FMC_time FMC_time_data



%% ---------------------------------------------------------------------- %
% Convert to dB and Plot                                                  %
% ---------------------------------------------------------------------- %%
tic;

if norm_to == 0
    max_ = 0;
    for view = 1 : Number_of_ims
        if max(abs(Sens(view).image(:))) > max_
            max_ = max(abs(Sens(view).image(:)));
        end
    end
else
    max_ = norm_to;
end

for view = 1 : Number_of_ims
   Sens(view).db_image = 20 * log10(abs(Sens(view).image) ./ max_); 
end

% Plot.
fig = figure(1);
% ax = repmat(subplot(plot_z, plot_x, 1), Number_of_ims, 1);
if image_block_info.type == "crack"
    sgtitle(sprintf('Sens %.2f Crack - %.2f deg', image_block_info.crack_length, rad2deg(image_block_info.angle)))
end
t = tiledlayout(plot_z, plot_x, 'TileSpacing', 'Compact');
for im = 1:Number_of_ims
    h(im) = nexttile;
%     ax(im) = subplot(plot_z, plot_x, im);
    imagesc(x*UC, z*UC, Sens(im).db_image);
    hold on
    title(Sens(im).name)
    caxis([-db_range_for_output, 0])
    plot(probe_coords(:, 1)*UC, probe_coords(:, 3)*UC, 'go');
    for wall = 1:size(geometry, 1)
        plot(geometry(wall).coords(:, 1)*UC, geometry(wall).coords(:, 3)*UC, 'r')
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



time_4 = double(toc);

fn_print_time('Plotted', time_4)

times = [time_1, time_2, time_3, time_4];



if savepath ~= ""
    cd(savepath)
    filename_fig = sprintf('%s.fig', savename);
    filename_mat = sprintf('%s.mat', savename);
    savefig(filename_fig)
    save(filename_mat, 'times', 'Sens', 'image_block_info')
end
% close all

% clear



end
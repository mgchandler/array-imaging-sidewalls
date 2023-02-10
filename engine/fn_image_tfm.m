function [Ims, Views_im, Views] = fn_image_tfm(FMC_time, FMC_time_data, model_options, varargin)

%% ---------------------------------------------------------------------- %
% Unpack model_config and model_options                                   %
% ---------------------------------------------------------------------- %%

PITCH = model_options.probe.width + model_options.probe.separation;
PIXEL = model_options.model.pixel;

probe_els = model_options.probe.num_els;
xmin = model_options.model.image_range(1);
xmax = model_options.model.image_range(2);
zmin = model_options.model.image_range(3);
zmax = model_options.model.image_range(4);
image_locs = model_options.model.image_locs;
scat_info = model_options.mesh.scat;
savepath = model_options.model.savepath;
savename = model_options.model.savename;
geometry = model_options.mesh.geom.geometry;
wall_for_imaging = model_options.model.wall_for_imaging;
boxsize = model_options.model.boxsize;
max_no_reflections = model_options.model.max_no_reflections;

probe_angle = model_options.probe.angle;
probe_standoff = model_options.probe.standoff;
probe_frequency = model_options.probe.freq;
el_length = model_options.probe.width;
couplant_speed = model_options.material.couplant_v;
couplant_density = model_options.material.couplant_density;
solid_long_speed = model_options.material.v_L;%sqrt(model_options.material.modulus * (1 - model_options.material.poisson) / (model_options.material.density * (1 + model_options.material.poisson) * (1 - 2*model_options.material.poisson)));
solid_shear_speed = model_options.material.v_S;%sqrt(model_options.material.modulus / (2 * model_options.material.density * (1 + model_options.material.poisson)));
solid_density = model_options.material.density;
max_num_reflections = model_options.model.max_no_reflections;
norm_to = model_options.model.norm_to;
db_range_for_output = model_options.model.db_range;
npw = model_options.mesh.n_per_wl;
time_it = model_options.model.time_it;
plot_it = model_options.model.plot_it;

if time_it
    tic
end

no_walls = size(geometry, 1);
no_cycles = model_options.probe.cycles;
frequency = model_options.probe.freq;

% Check whether we are in contact or immersion. If we are in contact, there
% will be no frontwall, and probe_standoff and probe_angle must equal zero.
% If we are in immersion, there must be a frontwall and probe_standoff must
% be non-zero.
is_frontwall = false;
for wall = 1:no_walls
    if geometry(wall).name == "F"
        is_frontwall = true;
        break
    end
end
if and(is_frontwall, probe_standoff ~= 0)
    is_contact = false;
elseif and(~is_frontwall, and(probe_standoff == 0, probe_angle == 0))
    is_contact = true;
else
    probe_standoff = .02e-3;
    is_contact = true;
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

num_paths = 0;
for num_reflections_in_path = 0:max_no_reflections
    num_paths = num_paths + 2^(num_reflections_in_path + 1);
end

Path_info_list = fn_direct_path_info(is_contact, probe_coords, model_options);
if max_no_reflections > 0
    Path_info_list = [Path_info_list; fn_skip_path_info(is_contact, probe_coords, model_options)];
end

%% ---------------------------------------------------------------------- %
% Scatterer Rays for box                                                  %
% ---------------------------------------------------------------------- %%

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
Views = fn_make_views(Paths, 1);

%% ---------------------------------------------------------------------- %
% Imaging Setup                                                           %
% ---------------------------------------------------------------------- %%
if time_it
    tic;
end

if image_locs == 0

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

else
    image_block = image_locs;
    im_x = unique(image_block(:, 1));
    im_z = unique(image_block(:, 3));
    xpts = size(im_x, 1)-1;
    zpts = size(im_z, 1)-1;
end

image_block_info = fn_scat_info("image", image_block(:, 1), image_block(:, 2), image_block(:, 3));
scatterer_coords = reshape(image_block, zpts+1, xpts+1, 3);
are_points_in_geometry = (scatterer_coords(:,:,1) >= xmin) .* ...
                         (scatterer_coords(:,:,1) <= xmax) .* ...
                         (scatterer_coords(:,:,3) >= zmin) .* ...
                         (scatterer_coords(:,:,3) <= zmax);
are_points_in_geometry = ones(size(are_points_in_geometry));
                         
% If Views_im passed in, do not recompute.
if nargin > 3
    for narg = 1 : nargin-3
        arg = varargin{narg};
        if isa(arg, 'struct')
            Views_im = arg;
        end
    end
else
    if time_it
        tic
    end

    Paths_im = repmat(fn_compute_ray(image_block_info, Path_info_list(1)), 1, num_paths);
    thispath = 2;
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
end

Number_of_ims = size(Views_im, 1);
% for im = 1:Number_of_ims
if isstruct(scat_info.fmc_mask)
    FMC_for_plotting = abs(FMC_time_data) .* scat_info.fmc_mask.data;
else
    FMC_for_plotting = abs(FMC_time_data);
end

% im_for_plotting = 1;
% fn_plot_FMC_at_time(FMC_for_plotting, FMC_time, Views_im(im_for_plotting).path_1.path_info, Views_im(im_for_plotting).path_2.path_info, [scat_info.x, scat_info.y, scat_info.z], sprintf('%s_(%s)_FMC.fig', savename, strrep(Views_im(im_for_plotting).name, ' ', '')));

Ims = repmat(fn_create_im("-", xpts+1, zpts+1), Number_of_ims, 1);
for view = 1 : Number_of_ims
    Ims(view).x = im_x;
    Ims(view).z = im_z;
    Ims(view).name = Views_im(view).name;
end

if time_it
    time_2 = double(toc);
    fn_print_time('Rays traced', time_2)
end

clear PIXEL xsize zsize num_paths Path_info_list image_block
clear image_block_info scatterer_coords Paths_im

% return



%% ---------------------------------------------------------------------- %
% Imaging                                                                 %
% ---------------------------------------------------------------------- %%
if time_it
    tic;
end

% Lookup times in FMC data.
if isstruct(scat_info.fmc_mask)
    FMC_time_data = FMC_time_data .* scat_info.fmc_mask.data;
end

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

if time_it
    time_4 = double(toc);
    fn_print_time('Imaged', time_4)
end

clear FMC_time FMC_time_data xpts zpts are_points_in_geometry tr_pair view tau



%% ---------------------------------------------------------------------- %
% Plotting.                                                               %
% ---------------------------------------------------------------------- %%
if time_it
    tic
end

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
%    Ims(view).phase_image = angle(Ims(view).image);
end

% im_idxs = [4, 12, 17];


% Plot.
for im = 1:Number_of_ims
    plot_idx = 1;
    Ims(im).plotExtras(plot_idx).x = probe_coords(:, 1);
    Ims(im).plotExtras(plot_idx).z = probe_coords(:, 3);
    Ims(im).plotExtras(plot_idx).color = 'g';
    Ims(im).plotExtras(plot_idx).marker = 'o';
    Ims(im).plotExtras(plot_idx).lineStyle = 'none';
    plot_idx = plot_idx + 1;
    for wall = 1:size(geometry, 1)
        % Only get first and last to reduce the amount saved - this will
        % need updating if we ever move away from polygonal geometry.
        Ims(im).plotExtras(plot_idx).x = geometry(wall).coords(:, 1);%[geometry(wall).coords(1, 1), geometry(wall).coords(end, 1)];
        Ims(im).plotExtras(plot_idx).z = geometry(wall).coords(:, 3);%[geometry(wall).coords(1, 3), geometry(wall).coords(end, 3)];
        Ims(im).plotExtras(plot_idx).color = 'r';
        Ims(im).plotExtras(plot_idx).marker = 'none';
        Ims(im).plotExtras(plot_idx).lineStyle = '-';
        plot_idx = plot_idx + 1;
    end
    
    % Plot scatterer
    Ims(im).plotExtras(plot_idx).x = scat_info.x;
    Ims(im).plotExtras(plot_idx).z = scat_info.z;
    Ims(im).plotExtras(plot_idx).color = 'r';
    Ims(im).plotExtras(plot_idx).marker = '.';
    Ims(im).plotExtras(plot_idx).lineStyle = 'none';
    plot_idx = plot_idx + 1;
    
    if boxsize ~= 0
        for s = 1 : size(scat_info.x, 1)
            if ~strcmp(scat_info.type, 'image')
%                 new_box_x = scat_info.x(s) + scat_info.r(s)/2*(sin(mean(Views(im).scat_inc_angles))+sin(mean(Views(im).scat_out_angles)));
%                 new_box_z = scat_info.z(s) + scat_info.r(s)/2*(cos(mean(Views(im).scat_inc_angles))+cos(mean(Views(im).scat_out_angles)));

                Ims(im).plotExtras(plot_idx).x = [scat_info.x(s) - boxsize/2, scat_info.x(s) - boxsize/2, scat_info.x(s) + boxsize/2, scat_info.x(s) + boxsize/2, scat_info.x(s) - boxsize/2];
                Ims(im).plotExtras(plot_idx).z = [scat_info.z(s) - boxsize/2, scat_info.z(s) + boxsize/2, scat_info.z(s) + boxsize/2, scat_info.z(s) - boxsize/2, scat_info.z(s) - boxsize/2];
                Ims(im).plotExtras(plot_idx).color = 'r';
                Ims(im).plotExtras(plot_idx).marker = 'none';
                Ims(im).plotExtras(plot_idx).lineStyle = '-';
                plot_idx = plot_idx + 1;
%                 Ims(im).plotExtras(plot_idx).x = [new_box_x - boxsize/2, new_box_x - boxsize/2, new_box_x + boxsize/2, new_box_x + boxsize/2, new_box_x - boxsize/2];
%                 Ims(im).plotExtras(plot_idx).z = [new_box_z - boxsize/2, new_box_z + boxsize/2, new_box_z + boxsize/2, new_box_z - boxsize/2, new_box_z - boxsize/2];
%                 Ims(im).plotExtras(plot_idx).color = 'g';
%                 Ims(im).plotExtras(plot_idx).marker = 'none';
%                 Ims(im).plotExtras(plot_idx).lineStyle = '-';
%                 plot_idx = plot_idx + 1;
                
                if isfield(Views(im).path_1, 'coords')
                    for leg = 1:size(Views(im).path_1.coords, 3)-1
                        Ims(im).plotExtras(plot_idx).x = [mean(Views(im).path_1.coords(:, s, leg, 1)), mean(Views(im).path_1.coords(:, s, leg+1, 1))];
                        Ims(im).plotExtras(plot_idx).z = [mean(Views(im).path_1.coords(:, s, leg, 3)), mean(Views(im).path_1.coords(:, s, leg+1, 3))];
                        Ims(im).plotExtras(plot_idx).color = [.5,.5,.5];
                        Ims(im).plotExtras(plot_idx).marker = 'none';
                        Ims(im).plotExtras(plot_idx).lineStyle = '-';
                        plot_idx = plot_idx + 1;
                    end
                end
                if isfield(Views(im).path_2, 'coords')
                    for leg = 1:size(Views(im).path_2.coords, 3)-1
                        Ims(im).plotExtras(plot_idx).x = [mean(Views(im).path_2.coords(:, s, leg, 1)), mean(Views(im).path_2.coords(:, s, leg+1, 1))];
                        Ims(im).plotExtras(plot_idx).z = [mean(Views(im).path_2.coords(:, s, leg, 3)), mean(Views(im).path_2.coords(:, s, leg+1, 3))];
                        Ims(im).plotExtras(plot_idx).color = [.5,.5,.5];
                        Ims(im).plotExtras(plot_idx).marker = 'none';
                        Ims(im).plotExtras(plot_idx).lineStyle = '-';
                        plot_idx = plot_idx + 1;
                    end
                end
            end
        end
    end
end

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

if plot_it
    fn_image_from_mat(Ims);
end

if time_it
    time_5 = double(toc);
    fn_print_time('Plotted', time_5)
end



if savepath ~= ""
%     cd(savepath)
    filename_fig = fullfile(savepath, sprintf('%s.fig', savename));
    filename_mat = fullfile(ssavepath, printf('%s.mat', savename));
    if plot_it
        savefig(filename_fig)
    end
    save(filename_mat, 'Ims')
end

end
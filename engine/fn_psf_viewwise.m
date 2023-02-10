function [Ims, Views_im, Views] = fn_psf_viewwise(model_options, varargin)
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

%% ---------------------------------------------------------------------- %
% Unpack model_config and model_options                                   %
% ---------------------------------------------------------------------- %%

savename = model_options.model.savename;
geometry = model_options.mesh.geom.geometry;

max_no_reflections = model_options.model.max_no_reflections;
pitch = model_options.probe.width + model_options.probe.separation;
probe_angle = model_options.probe.angle;
probe_els = model_options.probe.num_els;
probe_standoff = model_options.probe.standoff;
no_cycles = model_options.probe.cycles;
frequency = model_options.probe.freq;
scat_info = model_options.mesh.scat;
time_it = model_options.model.time_it;

model_options.model.time_it = 0;

if time_it
    tic
end

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
    is_contact = false;
elseif and(~is_frontwall, and(probe_standoff == 0, probe_angle == 0))
    is_contact = true;
else
    error('fn_sens: Invalid setup.')
end

rot_matrix = [cos(probe_angle) 0 sin(probe_angle); 0 1 0; -sin(probe_angle) 0 cos(probe_angle)];

probe_coords = zeros(3, probe_els);
probe_coords(1, :) = linspace(0, (probe_els - 1) * pitch, probe_els);
probe_coords(1, :) = probe_coords(1, :) - mean(probe_coords(1, :));
probe_coords = probe_coords.' * rot_matrix;
probe_coords(:, 3) = probe_coords(:, 3) - probe_standoff;

%% Precompute views
num_paths = 0;
for num_reflections_in_path = 0:max_no_reflections
    num_paths = num_paths + 2^(num_reflections_in_path + 1);
end
Path_info_list = fn_direct_path_info(is_contact, probe_coords, model_options);
if model_options.model.max_no_reflections > 0
    Path_info_list = [Path_info_list; fn_skip_path_info(is_contact, probe_coords, model_options)];
end

wall_for_imaging = model_options.model.wall_for_imaging;
frequency = model_options.probe.freq;
Paths = repmat(fn_compute_ray(model_options.mesh.scat, Path_info_list(1), geometry, frequency), 1, num_paths);
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
Views = fn_make_views(Paths, 1);

number_of_ims = size(Views, 1);

for im = 1:number_of_ims
    % Pass in the relevant view and paths only.
    model_options.model.model_views = Views(im);
    path1 = Views(im).path_1.path_info;
    path2 = Views(im).path_2.path_info;
    if strcmp(path1.name, path2.name)
        Path_info_list = [path1];
    else
        Path_info_list = [path1, path2];
    end
    model_options.model.path_info_list = Path_info_list;
    
    % Compute the FMC and TFM. Should produce only one image.
    [FMC_time, FMC_time_data] = fn_simulate_fmc(model_options);
    if im == 1
        [Im_view, Views_im, ~] = fn_image_tfm(FMC_time, FMC_time_data, model_options);
        Im_view = Im_view(1);
    else
        [Im_view, ~, ~] = fn_image_tfm(FMC_time, FMC_time_data, model_options, Views_im(im));
    end
    close all
    
    % Collect all into one.
    if im == 1
        Ims = repmat(Im_view, number_of_ims, 1);
    else
        Ims(im) = Im_view;
    end
end
if time_it
    time_2 = double(toc);
    fn_print_time('PSF calculated', time_2)
end

% fn_image_from_mat(Ims)

end
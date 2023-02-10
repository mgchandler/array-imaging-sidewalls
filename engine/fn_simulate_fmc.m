function [FMC_time, FMC_time_data] = fn_simulate_fmc(model_options)
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

PITCH = model_options.probe.width + model_options.probe.separation;

probe_els = model_options.probe.num_els;
scat_info = model_options.mesh.scat;
savepath = model_options.model.savepath;
savename = model_options.model.savename;
geometry = model_options.mesh.geom.geometry;
wall_for_imaging = model_options.model.wall_for_imaging;

probe_angle = model_options.probe.angle;
probe_standoff = model_options.probe.standoff;
probe_frequency = model_options.probe.freq;
el_length = model_options.probe.width;
couplant_speed = model_options.material.couplant_v;
couplant_density = model_options.material.couplant_density;
solid_long_speed = model_options.material.v_L;
solid_shear_speed = model_options.material.v_S;
solid_density = model_options.material.density;
max_no_reflections = model_options.model.max_no_reflections;
max_no_reverberations = model_options.model.max_no_reverberations;
model_geometry = model_options.model.model_geom;
multi_freq = model_options.model.multi_freq;
npw = model_options.mesh.n_per_wl;
time_it = model_options.model.time_it;

if time_it
    tic
end

% Additional parameters not directly dependent on inputs.
oversampling = 4;
no_cycles = 5;

no_walls = size(geometry, 1);

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
%     error('fn_sens: Invalid setup.')
end

clear wall



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



%% ---------------------------------------------------------------------- %
% Compute all path info                                                   %
% ---------------------------------------------------------------------- %%

if ~isstruct(model_options.model.model_views)
    num_paths = 0;
    for num_reflections_in_path = 0:max_no_reflections
        num_paths = num_paths + 2^(num_reflections_in_path + 1);
    end
    
    Path_info_list = fn_direct_path_info(is_contact, probe_coords, model_options);
    if max_no_reflections > 0
        Path_info_list = [Path_info_list; fn_skip_path_info(is_contact, probe_coords, model_options)];
    end



    %% ---------------------------------------------------------------------- %
    % Create simulation views                                                 %
    % ---------------------------------------------------------------------- %%
    % Precompute views for all scatterers, as well as imaging views.
    if time_it
        tic;
    end
    
    % Compute Imaging Paths
    if multi_freq
    %     error("fn_simulate_fmc: multi frequency not implemented")
        Paths = repmat(fn_compute_ray(scat_info, Path_info_list(1), geometry), 1, num_paths);
        path = 1;
        ii = 1;
        while path < num_paths
            ii = ii + 1;
            if length(Path_info_list(ii).speeds) > 1
                if wall_for_imaging == Path_info_list(ii).path_geometry.name
                    path = path + 1;
                    Paths(path) = fn_compute_ray(scat_info, Path_info_list(ii), geometry);
                end
            else % Must be direct
                path = path+1;
                Paths(path) = fn_compute_ray(scat_info, Path_info_list(ii), geometry);
            end
        end
    else
        frequency = probe_frequency;
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
    end
    
    % Create views from these paths.
    Views = fn_make_views(Paths, 0);
else
    Views = model_options.model.model_views;
end

if time_it
    time_1 = double(toc);
    fn_print_time('Simulation setup', time_1)
end



%% ---------------------------------------------------------------------- %
% Geometry setup.   NEEDS WORK                                            %
% ---------------------------------------------------------------------- %%

if model_geometry

    backwall_views = fn_make_geometry_views(probe_coords, geometry, ...
        [couplant_speed solid_long_speed solid_shear_speed], ...
        [couplant_density solid_density], probe_frequency, PITCH, el_length, ...
        max_no_reverberations, npw, is_contact, is_frontwall ...
    );
end

if time_it
    time_2 = double(toc);
    fn_print_time('Rays traced', time_2)
end

clear probe_angle el_length couplant_density solid_density PITCH



%% ------------------------------------------------------------------------
% Input signal
% -------------------------------------------------------------------------

% Create input signal.
time_step = 1 / (probe_frequency * oversampling);
max_t = 0;
for view = 1:size(Views, 1)
    max_in_view = max(Views(view).min_times(:));
    max_t = max(max_t, max_in_view);
end
if model_geometry
    for view = 1:size(backwall_views, 1)
        max_in_view = max(backwall_views(view).min_times(:));
        max_t = max(max_t, max_in_view);
    end
end
max_t = 1.25 * max_t; % Bit of extra room on the end
time_pts = ceil(max_t / time_step)+1;
[t11, t12, freq, in_freq_spec, fft_pts] = fn_create_input_signal(time_pts, probe_frequency, time_step , no_cycles);

clear probe_standoff oversampling time_step



%% ---------------------------------------------------------------------- %
% Multi-frequency ray weights                                             %
% ---------------------------------------------------------------------- %%

if multi_freq
    if time_it
        tic
    end
    frequency = freq(freq ~= 0);
    scat_info = fn_scat_info( ...
        model_options.mesh.scat.type, ...
        model_options.mesh.scat.x, ...
        model_options.mesh.scat.y, ...
        model_options.mesh.scat.z, ...
        model_options.mesh.scat.r, ...
        solid_long_speed  ./ frequency, ...
        solid_shear_speed ./ frequency, ...
        model_options.mesh.scat.angle, ...
        'ang_pts_over_2pi', 120 ...
    );
    for path = 1:num_paths
        Paths(path).scat_info = scat_info;
        Paths(path).weights = fn_compute_ray_weights(Paths(path), frequency, geometry);
        Paths(path).freq_array = frequency;
    end
    % Create views from these paths.
    Views = fn_make_views(Paths, 0);
    if time_it
        time_2b = double(toc);
        fn_print_time('Model coeffs computed', time_2b)
    end
end



%% ---------------------------------------------------------------------- %
% Simulation.                                                             %
% ---------------------------------------------------------------------- %%
if time_it
    tic
end

% Scatterer simulation.
out_freq_spec = 0;
num_scatterers = size(scat_info.x, 1);

for scatterer = 1 : num_scatterers
    for view = 1 : size(Views, 1)
        weights = Views(view).weights(:, scatterer, :);
        scat_amp = Views(view).scat_amps(:, scatterer, :);
        valid_path = Views(view).valid_path(:, scatterer);
        amp = (scat_amp .* weights);% .* valid_path);
        out_freq_spec = out_freq_spec + ...
            fn_propagate_spectrum(freq, in_freq_spec, Views(view).min_times(:, scatterer), amp, 0);

        clear weights scat_amp amp
    end
end

if model_geometry
    [num_geom_views, ~] = size(backwall_views);
    for view = 1 : num_geom_views
        bw_amp = conj(backwall_views(view).weights);% .* backwall_views(view).valid_path);
        out_freq_spec = out_freq_spec + ...
            fn_propagate_spectrum(freq, in_freq_spec, backwall_views(view).min_times, bw_amp, 0);
    end
end

% Convert back to time.
[FMC_time, FMC_time_data] = fn_convert_spectrum_to_time(freq, out_freq_spec, fft_pts, time_pts);
FMC_time = FMC_time';

% Hilbert Filtering (?) from fast_DAS function
diagonals = spdiags([1:length(FMC_time)]' < length(FMC_time)/2, 0, length(FMC_time), length(FMC_time));
FMC_time_data = ifft(diagonals * fft(FMC_time_data));

%     figure(3)
%     fn_plot_FMC_at_time(FMC_time_data, FMC_time, Path_info_list(3), Path_info_list(5), [[0, 0, 12.5e-3]], sprintf('%s_FMC.png', savename));

if time_it
    time_3 = double(toc);
    fn_print_time('Simulated', time_3)
end

if savepath ~= ""
    cd(savepath)
    filename =  sprintf('%s_FMC.mat', savename);
    save(filename, "FMC_time", "FMC_time_data", "Views")
end

end
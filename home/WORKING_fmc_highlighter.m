dirname = "C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Experimental Data\2022\L-shape FMCs with 1mm SDHs\22.09.12 - test + ind sample\FE";
chdir(dirname)

yaml_name = "a_0mm_d_1.yml";

fmc_name  = "a_0mm_d_1_1_BP.mat";

yaml_options = yaml.loadFile(yaml_name);

v_L = 6334.87;
v_S = 3113.59;

load(fmc_name);

yaml_options.data.time = time(:, 1);
yaml_options.data.data = data;

yaml_options.model.wall_for_imaging = "S1";
yaml_options.model.number_of_reflections = 1;

for kk = 1:size(yaml_options.mesh.geom.z, 2)
    yaml_options.mesh.geom.z{kk} = -yaml_options.mesh.geom.z{kk};
end
for kk = 1:size(yaml_options.mesh.scat.z, 2)
    yaml_options.mesh.scat.z{kk} = -yaml_options.mesh.scat.z{kk};
end
yaml_options.mesh.n_per_wl = 0;
yaml_options.model.pixel = 2.0e-3;

model_options = fn_default_model_options(yaml_options);

couplant_speed = model_options.material.couplant_v;
solid_long_speed = model_options.material.v_L;
solid_shear_speed = model_options.material.v_S;
couplant_density = model_options.material.couplant_density;
solid_density = model_options.material.density;
npw = model_options.mesh.n_per_wl;
probe_frequency = model_options.probe.freq;
PITCH = model_options.probe.width + model_options.probe.separation;
el_length = model_options.probe.width;
geometry = model_options.mesh.geom.geometry;
max_num_reflections = model_options.model.number_of_reflections;
PIXEL = model_options.model.pixel;

rot_matrix = [cos(model_options.probe.angle) 0 sin(model_options.probe.angle); 0 1 0; -sin(model_options.probe.angle) 0 cos(model_options.probe.angle)];
probe_coords = zeros(3, model_options.probe.num_els);
probe_coords(1, :) = linspace(0, (model_options.probe.num_els - 1) * model_options.probe.separation + model_options.probe.width, model_options.probe.num_els);
probe_coords(1, :) = probe_coords(1, :) - mean(probe_coords(1, :));
probe_coords = probe_coords.' * rot_matrix;
probe_coords(:, 3) = probe_coords(:, 3) - model_options.probe.standoff;
no_walls = size(geometry, 1);

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
        if strcmp(geometry(wall).name, model_options.model.wall_for_imaging)
            for mode1 = 0:1 % Mode of the first leg
                mode1_name = mode_names(mode1+1);
                for mode2 = 0:1 % Mode of the second leg
                    mode2_name = mode_names(mode2+1);
                    Skip_path_info = fn_path_info( ...
    ...%                %%  Use these names to include the wall in the path name
                            sprintf("%s %s %s", mode1_name, path_geometry.name, mode2_name), ...
                            sprintf("%s %s %s", mode2_name, path_geometry.name, mode1_name), ...
    ...%                %%  Use these names to exclude the wall in the path name
    ...%                     sprintf("%s %s", mode1_name, mode2_name), ...
    ...%                     sprintf("%s %s", mode2_name, mode1_name), ...
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
end

x_range = 30e-3:PIXEL:48e-3;
z_range = 26e-3:PIXEL:45e-3;
[x_range, z_range] = meshgrid(x_range, z_range);
image_coords = [x_range(:), zeros(length(x_range(:)), 1), z_range(:)];

savename = "NLoS region all full skip views"; % strcat("NLoS region ", tx_path.name, " - ", rx_path.rev_name, ".fig");
fn_plot_FMC_at_time(data, time(:, 1), Path_info_list(3:6), Path_info_list(3:6), image_coords, savename);

% for tx = 3:6
%     for rx = 3:6
%         tx_path = Path_info_list(tx);
%         rx_path = Path_info_list(rx);
%         savename = ""; % strcat("NLoS region ", tx_path.name, " - ", rx_path.rev_name, ".fig");
%         fn_plot_FMC_at_time(data, time(:, 1), tx_path, rx_path, image_coords, savename);
%     end
% end
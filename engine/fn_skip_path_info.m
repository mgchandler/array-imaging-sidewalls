function Path_info_list = fn_skip_path_info(is_contact, probe_coords, model_options)

mode_names = ["L", "T"];
speeds = [model_options.material.v_L, model_options.material.v_S];

num_paths = 4;
geometry = model_options.mesh.geom.geometry;
no_walls = size(geometry, 1);

Path_info_list = repmat(fn_path_info( ...
        "", ...
        "", ...
        [0], ...
        0, ...
        0, ...
        [model_options.material.couplant_v, model_options.material.v_L, model_options.material.v_S], ...
        [1], ...
        [model_options.material.couplant_density, model_options.material.density], ...
        0, ...
        0, ...
        0, ...
        probe_coords, ...
        model_options.mesh.n_per_wl ...
    ), ...
    num_paths, 1 ...
);

path = 1;
%% Contact
if is_contact
    for wall = 1:no_walls
        path_geometry = geometry(wall);
        for mode1 = 0:1 % Mode of the first leg
            mode1_name = mode_names(mode1+1);
            for mode2 = 0:1 % Mode of the second leg
                mode2_name = mode_names(mode2+1);
                wall_name = char(path_geometry.name);
                Skip_path_info = fn_path_info( ...
...%                %%  Use these names to include the wall in the path name
                    sprintf("%s %s %s", mode1_name, wall_name, mode2_name), ...
                    sprintf("%s %s %s", mode2_name, wall_name, mode1_name), ...
...%                %%  Use these names to exclude the wall in the path name
...%                         sprintf("%s %s", mode1_name, mode2_name), ...
...%                         sprintf("%s %s", mode2_name, mode1_name), ...
                    [mode1, mode2], ...
                    path_geometry, ...
                    [speeds(mode1+1), speeds(mode2+1)], ...
                    [model_options.material.couplant_v, model_options.material.v_L, model_options.material.v_S], ...
                    [1, 1], ...
                    [model_options.material.couplant_density, model_options.material.density], ...
                    model_options.probe.freq, ...
                    model_options.probe.width + model_options.probe.separation, ...
                    model_options.probe.width, ...
                    probe_coords, ...
                    model_options.mesh.n_per_wl ...
                );
                Path_info_list(path) = Skip_path_info;
                path = path + 1;
                clear Skip_path_info
            end
        end

        clear path_geometry

    end
%% Immersion
elseif ~is_contact
    wall_names = repmat("", size(geometry, 1), 1);
    for wall = 1:size(geometry, 1)
        wall_names(wall) = geometry(wall).name;
    end
    where_F = logical(count(wall_names, "F"));

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
                wall_name = chr(path_geometry(2).name);
                Skip_path_info = fn_path_info( ...
...%                %%  Use these names to include the wall in the path name
                    sprintf("%s %s %s", mode1_name, wall_name(1), mode2_name), ...
                    sprintf("%s %s %s", mode2_name, wall_name(1), mode1_name), ...
...%                %%  Use these names to exclude the wall in the path name
...%                         sprintf("%s %s", mode1_name, mode2_name), ...
...%                         sprintf("%s %s", mode2_name, mode1_name), ...
                    [0, mode1, mode2], ...
                    path_geometry, ...
                    [model_options.material.couplant_v, speeds(mode1+1), speeds(mode2+1)], ...
                    [model_options.material.couplant_v, model_options.material.v_L, model_options.material.v_S], ...
                    [0, 1, 1], ...
                    [model_options.material.couplant_density, model_options.material.density], ...
                    model_options.probe.freq, ...
                    model_options.probe.width + model_options.probe.separation, ...
                    model_options.probe.width, ...
                    probe_coords, ...
                    model_options.mesh.n_per_wl ...
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
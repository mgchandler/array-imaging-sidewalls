close all
clear
clc

dirname = "C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Experimental Data\2022\Amplitude Distribution Study";
cd(dirname)

yaml_name = "exp_scat_dist.yml";
yaml_options = yaml.loadFile(yaml_name);

b_or_s = "big";
suffix = " ACLs";

N = 1000;
psf_step = 5e-3;
boxsize_2 = 4.5e-3;

% FE
v_L = 6334.87;
v_S = 3113.59;

for kk = 1:size(yaml_options.mesh.geom.z, 2)
    yaml_options.mesh.geom.z{kk} = -yaml_options.mesh.geom.z{kk};
end
for kk = 1:size(yaml_options.mesh.scat.z, 2)
    yaml_options.mesh.scat.z{kk} = -yaml_options.mesh.scat.z{kk};
end
% yaml_options.mesh.scat = fn_scat_info( ...
%     "point", ...
%     0, 0, 0, ...
%     v_L, v_S ...
% );
yaml_options.model.time_it = false;
yaml_options.model.pixel = .5e-3;
yaml_options.mesh.geom.n_pts = 2000;
yaml_options.mesh.n_per_wl = 0;
yaml_options.model.max_no_reflections = 1;
yaml_options.model.savepath = "";%"C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\array-imaging-sidewalls\home-output\scat vs artefact distribution\RT Simulation\Random";
yaml_options.model.interp_method = "default";
yaml_options.model.time_it = 0;
yaml_options.model.model_geom = 0;
yaml_options.model.plot_it = 0;

%% Run through data sets
folders = dir(dirname);
ff = 1;
dd = 1;
% Sample m data sets randomly.
m = 10;
t1 = tic;
rng(mod(t1, 2^32))
samples = sort(floor(200*random('Uniform', 0, 1, m, 1)))+100;

stat_folders = ["non-defective - jig at corner - 0mm offset", "non-defective - jig at corner - 1mm offset", "non-defective - jig at corner 2 - 0mm offset"];
for ii = 1:length(folders)
    if and(and(folders(ii).isdir, and(~strcmp(folders(ii).name, '.'), ~strcmp(folders(ii).name, '..'))), any(strcmp(folders(ii).name, stat_folders)))
        thisdir = fullfile(dirname, folders(ii).name, 'TFMs', 'Relative coords');
        yaml_options.model.savepath = "";%thisdir;
        yaml_options.model.savename = "";
        for jj = 0:99
            try
                if ~any(jj == samples)
                    continue
                end
                cd(thisdir)
                matname = sprintf("%02d", jj);
%                 yaml_options.model.savename = sprintf("%s %s", matname, b_or_s);
                load(sprintf("%s sml sizing", matname))
                load(sprintf("%s %s", matname, b_or_s))
                % Update geometry. Do this manually as there's a
                % combination of mins/maxes which isn't easy to automate.
                yaml_options.mesh.geom.x{1} = max(geometry(1).point1(1), geometry(1).point2(1));
                yaml_options.mesh.geom.z{1} = min(geometry(1).point1(3), geometry(1).point2(3));
                yaml_options.mesh.geom.x{2} = max(geometry(2).point1(1), geometry(2).point2(1));
                yaml_options.mesh.geom.z{2} = max(geometry(2).point1(3), geometry(2).point2(3));
                yaml_options.mesh.geom.x{3} = max(geometry(3).point1(1), geometry(3).point2(1));
                yaml_options.mesh.geom.z{3} = max(geometry(3).point1(3), geometry(3).point2(3));
                yaml_options.mesh.geom.x{4} = max(geometry(4).point1(1), geometry(4).point2(1));
                yaml_options.mesh.geom.z{4} = max(geometry(4).point1(3), geometry(4).point2(3));
                yaml_options.mesh.geom.x{5} = max(geometry(5).point1(1), geometry(5).point2(1));
                yaml_options.mesh.geom.z{5} = max(geometry(5).point1(3), geometry(5).point2(3));
                yaml_options.mesh.geom.x{6} = max(geometry(5).point1(1), geometry(5).point2(1));
                yaml_options.mesh.geom.z{6} = min(geometry(5).point1(3), geometry(5).point2(3));
                yaml_options.model.image_range = [geometry(3).point1(1)+1e-5, geometry(1).point1(1)-1e-5, ...
                                                  geometry(4).point1(3)+1e-5, geometry(2).point1(3)-1e-5];
                model_options = fn_default_model_options(yaml_options);
                %% Step 1: Calculate PSF at range of locations. Do this only in lower bit of geometry.
                % Coordinates at which we find hs. Quite low res as it
                % varies slowly.

                psf_x = linspace(geometry(3).point1(1)+.01e-3, geometry(1).point1(1)-.01e-3, 8);%geometry(3).point1(1):psf_step:geometry(1).point1(1);
                psf_z = linspace(geometry(4).point1(3)+.01e-3, geometry(2).point1(3)-.01e-3, 8);%geometry(4).point1(3):psf_step:geometry(2).point1(3);
                hs = zeros(size(psf_x, 2), size(psf_z, 2), 21);
                t1 = tic;
                for x = 1:size(psf_x, 2)
                    for z = 1:size(psf_z, 2)
                        % Only do a small image to save runtime.
                        model_options.mesh.scat.x = psf_x(x);
                        model_options.mesh.scat.z = psf_z(z);
                        model_options.model.image_range = [psf_x(x) - boxsize_2, psf_x(x) + boxsize_2, psf_z(z) - boxsize_2, psf_z(z) + boxsize_2];
                        % Find it viewwise.
                        [Psf, ~, ~] = fn_psf_viewwise(model_options);
                        % Compute hs for each view at this location.
                        for im = 1:21
                            hs(x, z, im) = fn_acl(Psf(im).x, Psf(im).z, abs(Psf(im).image));
                        end
                    end
                end
                %% Step 2: Calculate h using hs found using psf.
                % We have a TFM of the full geometry, but we only really
                % care to image in the lower bit. Get the subset of the
                % TFM.
                [~, x_min] = min(abs(Ims(1).x - geometry(3).point1(1)), [], 'all');
                [~, x_max] = min(abs(Ims(1).x - geometry(1).point1(1)), [], 'all');
                [~, z_min] = min(abs(Ims(1).z - geometry(4).point1(3)), [], 'all');
                [~, z_max] = min(abs(Ims(1).z - geometry(2).point1(3)), [], 'all');
                x_min = x_min + 1;
                x_max = x_max - 1;
                z_min = z_min + 1;
                z_max = z_max - 1;
                x_step = 2*(Ims(1).x(2) - Ims(1).x(1));
                z_step = 2*(Ims(1).z(2) - Ims(1).z(1));
                h = zeros(x_max-x_min+1, z_max-z_min+1, 21);
                for x = x_min:x_max
                    for z = z_min:z_max
                        % We only compute ACF in a box-windowed region,
                        % whose size comes from hs. Find out what this is
                        % in array indices.
                        hs_box = interp2(psf_x, psf_z, squeeze(hs(:, :, im)), Ims(im).x(x), Ims(im).z(z), 'cubic');
                        x_range = max(x - ceil(hs_box/x_step), 1):min(x + ceil(hs_box/x_step), size(Ims(im).x, 2));
                        z_range = max(z - ceil(hs_box/z_step), 1):min(z + ceil(hs_box/z_step), size(Ims(im).z, 2));
                        for im = 1:21
                            h(x-x_min+1, z-z_min+1, im) = fn_acl(Ims(im).x(x_range), Ims(im).z(z_range), abs(Ims(im).image(z_range, x_range)));
                        end
                    end
                end
                save(sprintf("%s%s", matname, suffix), "h", "hs")
                time_2 = double(toc(t1));
                fn_print_time(sprintf('Mask %s calculated', matname), time_2)
            catch ME
                disp(ME.message)
            end
        end
        samples = samples - 100;
    end
end
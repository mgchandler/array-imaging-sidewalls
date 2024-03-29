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

percentile_target = .99;

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
samples = sort(floor(200*random('Uniform', 0, 1, m, 1)));

Masks = repmat(struct("matname", 0, "x", 0, "z", 0, "mask", 0, "geometry", 0), m, 1);
idx = 1;

min_x = 0;
min_z = 0;

figure
hold on

stat_folders = ["non-defective - jig at corner - 0mm offset", "non-defective - jig at corner - 1mm offset", "non-defective - jig at corner 2 - 0mm offset"];
for ii = 1:length(folders)
    if and(and(folders(ii).isdir, and(~strcmp(folders(ii).name, '.'), ~strcmp(folders(ii).name, '..'))), any(strcmp(folders(ii).name, stat_folders)))
        thisdir = fullfile(dirname, folders(ii).name, 'TFMs', 'Relative coords');
        yaml_options.model.savepath = "";%thisdir;
        yaml_options.model.savename = "";
        for jj = 0:99
            try
                %% Generate ACLs for ACF Mask
%                 if ~any(jj == samples)
%                     continue
%                 end
                cd(thisdir)
                matname = sprintf("%02d", jj);
% %                 yaml_options.model.savename = sprintf("%s %s", matname, b_or_s);
                load(sprintf("%s sml sizing", matname))
                load(sprintf("%s %s", matname, b_or_s))
%                 % Update geometry. Do this manually as there's a
%                 % combination of mins/maxes which isn't easy to automate.
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
%                 model_options = fn_default_model_options(yaml_options);
%                 %% Step 1: Calculate PSF at range of locations. Do this only in lower bit of geometry.
%                 % Coordinates at which we find hs. Quite low res as it
%                 % varies slowly.
% 
%                 psf_x = linspace(geometry(3).point1(1)+.01e-3, geometry(1).point1(1)-.01e-3, 8);%geometry(3).point1(1):psf_step:geometry(1).point1(1);
%                 psf_z = linspace(geometry(4).point1(3)+.01e-3, geometry(2).point1(3)-.01e-3, 8);%geometry(4).point1(3):psf_step:geometry(2).point1(3);
%                 hs = zeros(size(psf_x, 2), size(psf_z, 2), 21);
%                 t1 = tic;
%                 for x = 1:size(psf_x, 2)
%                     for z = 1:size(psf_z, 2)
%                         % Only do a small image to save runtime.
%                         model_options.mesh.scat.x = psf_x(x);
%                         model_options.mesh.scat.z = psf_z(z);
%                         model_options.model.image_range = [psf_x(x) - boxsize_2, psf_x(x) + boxsize_2, psf_z(z) - boxsize_2, psf_z(z) + boxsize_2];
%                         % Find it viewwise.
%                         [Psf, ~, ~] = fn_psf_viewwise(model_options);
%                         % Compute hs for each view at this location.
%                         for im = 1:21
%                             hs(x, z, im) = fn_acl(Psf(im).x, Psf(im).z, abs(Psf(im).image));
%                         end
%                     end
%                 end
%                 %% Step 2: Calculate h using hs found using psf.
%                 % We have a TFM of the full geometry, but we only really
%                 % care to image in the lower bit. Get the subset of the
%                 % TFM.
                [~, x_min] = min(abs(Ims(1).x - geometry(3).point1(1)), [], 'all');
                [~, x_max] = min(abs(Ims(1).x - geometry(1).point1(1)), [], 'all');
                [~, z_min] = min(abs(Ims(1).z - geometry(4).point1(3)), [], 'all');
                [~, z_max] = min(abs(Ims(1).z - geometry(2).point1(3)), [], 'all');
                x_min = x_min + 1;
                x_max = x_max - 1;
                z_min = z_min + 1;
                z_max = z_max - 1;
%                 x_step = 2*(Ims(1).x(2) - Ims(1).x(1));
%                 z_step = 2*(Ims(1).z(2) - Ims(1).z(1));
%                 h = zeros(x_max-x_min+1, z_max-z_min+1, 21);
%                 acl_x = Ims(1).x(x_min:x_max);
%                 acl_z = Ims(1).z(z_min:z_max);
%                 for x = x_min:x_max
%                     for z = z_min:z_max
%                         % We only compute ACF in a box-windowed region,
%                         % whose size comes from hs. Find out what this is
%                         % in array indices.
%                         hs_box = interp2(psf_x, psf_z, squeeze(hs(:, :, im)), Ims(im).x(x), Ims(im).z(z), 'cubic');
%                         x_range = max(x - ceil(hs_box/x_step), 1):min(x + ceil(hs_box/x_step), size(Ims(im).x, 2));
%                         z_range = max(z - ceil(hs_box/z_step), 1):min(z + ceil(hs_box/z_step), size(Ims(im).z, 2));
%                         for im = 1:21
%                             h(x-x_min+1, z-z_min+1, im) = fn_acl(Ims(im).x(x_range), Ims(im).z(z_range), abs(Ims(im).image(z_range, x_range)));
%                         end
%                     end
%                 end
%                 save(sprintf("%s%s", matname, suffix), "h", "hs", "psf_x", "psf_z", "acl_x", "acl_z")
%                 time_2 = double(toc(t1));
%                 fn_print_time(sprintf('Mask %s calculated', matname), time_2)
            %% Make ACF mask
%             load(sprintf("%s%s", matname, suffix))
%             [X, Z] = meshgrid(acl_x, acl_z);
%             h = permute(h, [2,1,3]);
%             hs_interp = zeros(size(h));
%             for im = 1:21
%                 hs_interp(:, :, im) = interp2(psf_z, psf_x, squeeze(hs(:, :, im)), Z, X);
%             end
%             acf_mask = h ./ hs_interp > 1;
%             %% Grow regions (HICS)
%             load(sprintf("%s %s", matname, b_or_s))
%             S = Ims;
%             [~, x_min] = min(abs(S(1).x - acl_x(1)));
%             [~, x_max] = min(abs(S(1).x - acl_x(end)));
%             [~, z_min] = min(abs(S(1).z - acl_z(1)));
%             [~, z_max] = min(abs(S(1).z - acl_z(end)));
%             mask = acf_mask;
%             for im = 1:21
%                 % Initialise S
%                 if im == 4
%                     a = 1;
%                 end
%                 S(im).x = acl_x;
%                 S(im).z = acl_z;
%                 S(im).image = S(im).image(z_min:z_max, x_min:x_max);
%                 S(im).db_image = 20 * log10(abs(S(im).image) ./ max(abs(S(im).image(:)), [], 'omitnan'));
%                 ratios = [];
%                 while true
%                     % Update S based on mask.
%                     S(im).image(squeeze(mask(:, :, im))) = nan;
%                     S(im).db_image(squeeze(mask(:, :, im))) = nan;
%                     % Fit plane and correct for it
%                     plane = fit([X(:), Z(:)], S(im).db_image(:), 'poly11', 'Exclude', reshape(mask(:, :, im), [], 1));
%                     c = plane.p10 * X + plane.p01 * Z;
%                     S(im).image = abs(S(im).image) ./ 10.^(c / 20);
%                     S(im).db_image = 20 * log10(abs(S(im).image) ./ max(abs(S(im).image(:)), [], 'omitnan'));
%                     % Fit Rayleigh to image
%                     i_c = 0:max(S(im).image(:), [], 'omitnan')/1000:max(S(im).image(:));
%                     try
%                         rayleigh = fitdist(S(im).image(:), 'Rayleigh');
%                     catch ME
%                         break
%                     end
%                     % Mark everything >99th percentile
%                     p = icdf(rayleigh, percentile_target);
%                     T = double(S(im).image > p);
%                     T(isnan(S(im).image)) = nan;
%                     ratio = sum(T(:), 'omitnan') / sum(~isnan(S(im).image(:)));
%                     if ratio < (1-percentile_target)
%                         break
%                     end
%                     % If there are NaNs in T, it is because they have already been masked.
%                     % Set to zero to aid plotting - mask already contains  1.
%                     T(isnan(T)) = 0;
%                     ratios = [ratios, ratio];
%                     % Pick the point which has the largest number of marked/masked neighbours, and mask that one.
%                     max_points = 0;
%                     max_loc = [0, 0];
%                     for x = 1:size(S(im).x, 2)
%                         for z = 1:size(S(im).z, 2)
%                             % If we have already masked this point, do not try to mask it again.
%                             if mask(z, x, im)
%                                 continue
%                             end
%                             % Find neighbourhood
%                             box_x = round(hs_interp(z, x, im) / (S(im).x(2) - S(im).x(1)));
%                             box_z = round(hs_interp(z, x, im) / (S(im).z(2) - S(im).z(1)));
%                             min_box_x = max(1, x-box_x);
%                             max_box_x = min(size(S(im).x, 2), x+box_x);
%                             min_box_z = max(1, z-box_z);
%                             max_box_z = min(size(S(im).z, 2), z+box_z);
%                             % Count neighbours
%                             neighbourhood = or(mask(min_box_z:max_box_z, min_box_x:max_box_x, im), T(min_box_z:max_box_z, min_box_x:max_box_x));
%                             if sum(neighbourhood(:)) / ((max_box_x - min_box_x)*(max_box_z - min_box_z)) > max_points
%                                 max_points = sum(neighbourhood(:)) / ((max_box_x - min_box_x)*(max_box_z - min_box_z));
%                                 max_loc = [x, z];
%                             end
%                         end
%                     end
% %                     figure(1)
% %                     title("T")
% %                     set(gcf, 'Position', [-3768, 355, 560, 420])
% %                     T(isnan(S(im).image)) = nan;
% %                     h = imagesc(T);
% %                     set(h, 'AlphaData', ~isnan(T))
% %                     figure(2)
% %                     title("~S")
% %                     set(gcf, 'Position', [-3192, 355, 560, 420])
% %                     imagesc(mask(:, :, im))
% %                     figure(3)
% %                     title("T / S")
% %                     set(gcf, 'Position', [-3404, -161, 560, 420])
% %                     hold off
% %                     plot(ratios)
% %                     hold on
% %                     plot([1, numel(ratios)], [1-percentile_target, 1-percentile_target])
% %                     figure(4)
% %                     title("Rayleigh Dist")
% %                     set(gcf, 'Position', [-2808, -161, 560, 420])
% %                     hold off
% %                     histogram(S(im).image)
% %                     hold on
% %                     plot(i_c, rayleigh.pdf(i_c))
% %                     plot([p, p], [0, max(rayleigh.pdf(i_c))])
% 
%                     mask(max_loc(2), max_loc(1), im) = true;
%                 end
%             end
%             save(sprintf("%s acf+hics mask", matname), "acf_mask", "mask", "acl_x", "acl_z")
            %% Simple thresholding mask
            load(sprintf("%s%s", matname, suffix))
            h = permute(h, [2,1,3]);
            [X, Z] = meshgrid(acl_x, acl_z);
            % Assume Rayleigh distribution
            S = Ims;
            [~, x_min] = min(abs(S(1).x - acl_x(1)));
            [~, x_max] = min(abs(S(1).x - acl_x(end)));
            [~, z_min] = min(abs(S(1).z - acl_z(1)));
            [~, z_max] = min(abs(S(1).z - acl_z(end)));
            mask = zeros(size(h));
            max_ = 0;
            for im = 1:21
                % Initialise S
                if im == 4
                    a = 1;
                end
                S(im).x = acl_x;
                S(im).z = acl_z;
                S(im).image = S(im).image(z_min:z_max, x_min:x_max);
                S(im).db_image = S(im).db_image(z_min:z_max, x_min:x_max);
%                 ratios = [];

                % Update S based on mask.
%                 S(im).image(squeeze(mask(:, :, im))) = nan;
%                 S(im).db_image(squeeze(mask(:, :, im))) = nan;
                % Fit plane and correct for it
%                 plane = fit([X(:), Z(:)], S(im).db_image(:), 'poly11', 'Exclude', reshape(mask(:, :, im), [], 1));
%                 c = plane.p10 * X + plane.p01 * Z;
%                 S(im).image = abs(S(im).image) ./ 10.^(c / 20);
%                 S(im).db_image = 20 * log10(abs(S(im).image) ./ max(abs(S(im).image(:)), [], 'omitnan'));
                % Fit Rayleigh to image
%                 i_c = 0:max(S(im).image(:), [], 'omitnan')/1000:max(S(im).image(:));
%                 try
%                 rayleigh = fitdist(S(im).image(:), 'Rayleigh');
%                 catch ME
%                     break
%                 end
                % Mark everything >99th percentile
%                 p = icdf(rayleigh, percentile_target);
                S(im).db_image = 20 * log10(abs(S(im).image) ./ max(abs(S(im).image(:)), [], 'omitnan'));
                T = double(S(im).db_image > -20);
                mask(:, :, im) = T;
            end

            save(sprintf("%s 20db_thresh mask", matname), "mask", "acl_x", "acl_z")

            %% Analyse all of the data generated.
% %             load(sprintf("%s acf+hics mask", matname))
%             disp(matname)
%             Ims(1).x = Ims(1).x(x_min:x_max);
%             Ims(1).z = Ims(1).z(z_min:z_max);
%             Ims(1).db_image = Ims(1).db_image(z_min:z_max, x_min:x_max);
% %             H   = repmat(struct("x", acl_x, "z", acl_z, "db_image", squeeze(h(:, :, 1)), "name", ""), 21, 1);
% %             Hs  = repmat(struct("x", acl_x, "z", acl_z, "db_image", squeeze(hs_interp(:, :, 1)), "name", ""), 21, 1);
% %             HHs = repmat(struct("x", acl_x, "z", acl_z, "db_image", squeeze(h(:, :, 1)) ./ squeeze(hs_interp(:, :, 1)), "name", ""), 21, 1);
% % 
% %             ACF  = repmat(struct("x", acl_x, "z", acl_z, "db_image", squeeze(acf_mask(:, :, 1)), "name", ""), 21, 1);
%             HICS = repmat(struct("x", acl_x, "z", acl_z, "db_image", squeeze(mask(:, :, 1)), "name", ""), 21, 1);
%             for im = 2:21
%                 Ims(im).x = Ims(im).x(x_min:x_max);
%                 Ims(im).z = Ims(im).z(z_min:z_max);
%                 Ims(im).db_image = Ims(im).db_image(z_min:z_max, x_min:x_max);
% % 
% %                 H(im).db_image   = squeeze(h(:, :, im));
% %                 Hs(im).db_image  = squeeze(hs_interp(:, :, im));
% %                 HHs(im).db_image = squeeze(h(:, :, im)) ./ squeeze(hs_interp(:, :, im));
% % 
% %                 ACF(im).db_image  = squeeze(acf_mask(:, :, im));
%                 HICS(im).db_image = squeeze(mask(:, :, im));
%             end
%             fn_image_from_mat(Ims)
%             grp = get(get(gcf, 'Children'), 'Children');
%             for im = 2:22
% %                 colormap gray
%                 grp(im).XLim = [1e3*Ims(1).x(1), 1e3*Ims(1).x(end)];
%                 grp(im).YLim = [1e3*Ims(1).z(1), 1e3*Ims(1).z(end)];
%             end
% 
% %             fn_image_from_mat(H)
% %             grp = get(get(gcf, 'Children'), 'Children');
% %             grp(1).Label.String = sprintf("%s h", matname);
% %             for im = 2:22
% % %                 colormap gray
% %                 grp(im).CLim = [0, .005];
% %             end
% %             fn_image_from_mat(Hs)
% %             grp = get(get(gcf, 'Children'), 'Children');
% %             grp(1).Label.String = sprintf("%s h_s", matname);
% %             for im = 2:22
% % %                 colormap gray
% %                 grp(im).CLim = [0, .008];
% %             end
% %             fn_image_from_mat(HHs)
% %             grp = get(get(gcf, 'Children'), 'Children');
% %             grp(1).Label.String = sprintf("%s h/h_s", matname);
% %             for im = 2:22
% % %                 colormap gray
% %                 grp(im).CLim = [0, 2];
% %             end
% %            
% %             fn_image_from_mat(ACF)
% %             grp = get(get(gcf, 'Children'), 'Children');
% %             grp(1).Label.String = sprintf("%s ACF mask", matname);
% %             for im = 2:22
% %                 colormap gray
% %                 grp(im).CLim = [0, 1];
% %             end
%             fn_image_from_mat(HICS)
%             grp = get(get(gcf, 'Children'), 'Children');
%             grp(1).Label.String = sprintf("%s ACF+HICS mask", matname);
%             for im = 2:22
%                 colormap gray
%                 grp(im).CLim = [0, 1];
%             end
            %% Combine all 10 masks into one image.
            mask_type = "20db_thresh mask";
            load(sprintf("%s %s", matname, mask_type))
            Masks(idx).geometry = geometry;
            Masks(idx).matname = matname;
            Masks(idx).x = acl_x - geometry(1).point1(1);
            Masks(idx).z = acl_z - geometry(4).point1(3);
            Masks(idx).mask = mask;

%             plot([geometry(1).point1(1), geometry(1).point1(1), geometry(3).point1(1), geometry(3).point1(1), geometry(5).point1(1)], ...
%                  [geometry(1).point1(3), geometry(2).point1(3), geometry(2).point1(3), geometry(4).point1(3), geometry(4).point1(3)], ...
%                  'r')

            idx = idx + 1;
            catch ME
                disp(ME.message)
            end

        end
%         samples = samples - 100;
    end
end

% Want to redefine acl_x and acl_z so that they are all on the same base.
% Subtract S2.point1.

% mask_type = "acf mask";
max_x = 0;
max_z = 0;
for mask = 1:m
    if size(Masks(mask).x, 2) > max_x
        max_x = size(Masks(mask).x, 2);
    end
    if size(Masks(mask).z, 2) > max_z
        max_z = size(Masks(mask).z, 2);
    end
end
global_mask = zeros(max_z, max_x, 21);
unmasked = zeros(max_z, max_x, 21);
for mask = 1:m
    global_mask(1:size(Masks(mask).z, 2), end-size(Masks(mask).x, 2)+1:end, :) = or(global_mask(1:size(Masks(mask).z, 2), end-size(Masks(mask).x, 2)+1:end, :), Masks(mask).mask);
    unmasked(1:size(Masks(mask).z, 2), end-size(Masks(mask).x, 2)+1:end, :) = or(unmasked(1:size(Masks(mask).z, 2), end-size(Masks(mask).x, 2)+1:end, :), ones(size(Masks(mask).z, 2), size(Masks(mask).x, 2), 21));
end
global_mask = and(global_mask, unmasked);
save(sprintf("%s global", mask_type), "global_mask")

Global = Ims;
for im = 1:21
    Global(im).db_image = global_mask(:, :, im);
    Global(im).x = acl_x;
    Global(im).z = acl_z;
end
Global = rmfield(Global, "plotExtras");
fn_image_from_mat(Global)
grp = get(get(gcf, 'Children'), 'Children');
for im = 2:22
    grp(im).CLim = [0, 1];
    colormap gray
end

% New_masks = Masks;
% for mask = 1:m
%     new_mask = zeros(max_z, max_x, 21);
%     new_x = linspace(Masks(mask).x(1), Masks(mask).x(end), max_x);
%     new_z = linspace(Masks(mask).z(1), Masks(mask).z(end), max_z);
%     [New_x, New_z] = meshgrid(new_x, new_z);
%     for im = 1:21
%         new_mask(:, :, im) = interp2(Masks(mask).x, Masks(mask).z, double(Masks(mask).mask(:, :, im)), New_x, New_z, 'cubic');
%     end
%     new_mask(new_mask > .5) = 1;
%     new_mask(new_mask <= .5) = 0;
%     New_masks(mask).x = new_x;
%     New_masks(mask).z = new_z;
%     New_masks(mask).mask = new_mask;
% end
% global_mask = zeros(max_z, max_x, 21);
% for mask = 1:m
%     this_mask = New_masks(mask).mask;
%     global_mask = or(global_mask, this_mask);
% end
% Global_mask = repmat(struct("x", 0, "z", 0, "db_image", 0), 21, 1);
% for im = 1:21
%     Global_mask(im).name = "";
%     Global_mask(im).x = 1:max_x;
%     Global_mask(im).z = 1:max_z;
%     Global_mask(im).db_image = global_mask(:, :, im);
% end
% save(sprintf("%s global", mask_type), "global_mask")
% fn_image_from_mat(Global_mask)
% grp = get(get(gcf, 'Children'), 'Children');
% grp(1).Label.String = "global mask";
% for im = 2:22
%     colormap gray
%     grp(im).CLim = [0, 1];
% end
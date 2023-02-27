close all
clear
clc

dirname = "C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Experimental Data\2022\Amplitude Distribution Study";
cd(dirname)

yaml_name = "exp_scat_dist.yml";
yaml_options = yaml.loadFile(yaml_name);

b_or_s = "big";
mask_type = "";

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

best_views = [];

%% Run through data sets
folders = dir(dirname);
ff = 1;
dd = 1;
% Sample m data sets randomly.
m = 10;
t1 = tic;
rng(mod(t1, 2^32))
samples = sort(floor(200*random('Uniform', 0, 1, m, 1)));

% load(sprintf("%s mask global", mask_type))
load("big rician params")
load(fullfile("defective - jig at corner - 0mm offset", "TFMs", "Relative coords", "01 big"))
load(fullfile("defective - jig at corner - 0mm offset", "TFMs", "Relative coords", "01 sml sizing"))
load(fullfile("defective - jig at corner - 0mm offset", "TFMs", "Relative coords", "01 sens"))
params_geom = geometry;
params_x = Ims(1).x - params_geom(1).point1(1);
params_z = Ims(1).z - params_geom(2).point1(3);
sigma = zeros(size(params_z, 2), size(params_x, 2), 21);
for im = 1:21
    sigma(:, :, im) = reshape(nd_rician(2, im, :), size(params_z, 2), size(params_x, 2));
    Sens(im).x = Sens(im).x - params_geom(1).point1(1);
    Sens(im).z = Sens(im).z - params_geom(2).point1(3);
end
idx = 1;
thresholds = linspace(0, .01, 10000);
% log10_thresholds = linspace(-300, 0, 1000);
% thresholds = [0, 10.^log10_thresholds];

max_i = 0;

stat_folders = ["non-defective - jig at corner - 0mm offset", "non-defective - jig at corner - 1mm offset", "non-defective - jig at corner 2 - 0mm offset"];
stat_folders = ["defective - jig at corner - 0mm offset", "defective - jig at corner - 1mm offset", "defective - jig at corner 2 - 0mm offset", "non-defective - jig at corner - 0mm offset", "non-defective - jig at corner - 1mm offset", "non-defective - jig at corner 2 - 0mm offset"];
for ii = 1:length(folders)
    if and(and(folders(ii).isdir, and(~strcmp(folders(ii).name, '.'), ~strcmp(folders(ii).name, '..'))), any(strcmp(folders(ii).name, stat_folders)))
        thisdir = fullfile(dirname, folders(ii).name, 'TFMs', 'Relative coords');
        yaml_options.model.savepath = "";%thisdir;
        yaml_options.model.savename = "";
        pass_fail = false(100, size(thresholds, 2));
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
%                 load(sprintf("%s %s", matname, b_or_s))
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
                model_options = fn_default_model_options(yaml_options);
                model_options.mesh.geom.geometry = geometry;
%                 [Ims, ~, ~] = fn_tfm(model_options);
%                 [Sens, ~] = fn_sens(model_options);
                %% Apply the mask. Scale it to the size of geometry.
%                 sx = geometry(3).point1(1);
%                 sz = geometry(3).point1(3);
%                 ex = geometry(1).point2(1);
%                 ez = geometry(1).point2(3);
%                 nx = size(global_mask, 2);
%                 nz = size(global_mask, 1);
%                 mask_x = linspace(sx, ex, nx);
%                 mask_z = linspace(sz, ez, nz);
%                 [im_x, im_z] = meshgrid(Ims(1).x, Ims(1).z);
%                 % Rescale the mask.
%                 scaled_mask = zeros(size(Ims(1).z, 2), size(Ims(1).x, 2), 21);
%                 for im = 1:21
%                     scaled_mask(:, :, im) = interp2(mask_x, mask_z, double(global_mask(:, :, im)), im_x, im_z);
%                 end
%                 scaled_mask(scaled_mask > .5) = 1;
%                 scaled_mask(scaled_mask <= .5) = 0;
                % Apply it.
%                 for im = 1:21
%                     Ims(im).image(logical(scaled_mask(:, :, im))) = nan;
%                     Ims(im).db_image(logical(scaled_mask(:, :, im))) = nan;
% %                     Ims(im).image = Ims(im).image .* scaled_mask(:, :, im);
% %                     Ims(im).db_image = Ims(im).db_image .* scaled_mask(:, :, im);
%                 end
                % Get the relevant bits of sigma.
                
                %% Do the ROC bit.
%                 % What is the best view? Within ROC, find where sensitivity is highest for pixel.
%                 % Find out where Sens and Ims overlap. 
%                 ims_x  = reshape(Ims(1).x, 1, []) - geometry(1).point1(1);
%                 ims_z  = reshape(Ims(1).z, 1, []) - geometry(2).point1(3);
%                 % We want the valid subset of each.
%                 [~, sx]  = min(abs(params_x - ims_x(1)));
%                 [~, sz]  = min(abs(params_z - ims_z(1)));
%                 [~, ex]  = min(abs(params_x - ims_x(end)));
%                 [~, ez]  = min(abs(params_z - ims_z(end)));
%                 sigma = zeros(ez-sz+1, ex-sx+1, 21);
%                 for im = 1:21
%                     Ims(im).x = ims_x;
%                     Ims(im).z = ims_z;
%                     Sens(im).x = ims_x;
%                     Sens(im).z = ims_z;
% 
%                     im_sigma = reshape(nd_rician(2, im, :), size(params_z, 2), size(params_x, 2));
%                     sigma(:, :, im) = im_sigma(sz:ez, sx:ex);
%                 end
%                 save(sprintf("%s large roc input", matname), "Sens", "Ims", "sigma")
%                 continue

                %% ROC analysis on masked data.
                load(sprintf("%s large roc input", matname))
                sx = min(size(Ims(1).x, 2));%, size(global_mask, 2));
                sz = min(size(Ims(1).z, 2));%, size(global_mask, 1));
% 
%                 for im = 1:21
%                     if max(abs(Ims(im).image(:))) > max_i
%                         max_i = max(abs(Ims(im).image(:)));
%                     end
%                 end

                % Mask each image, and work out best view.
                metric_best = zeros(sz, sx);
                metric_arg = zeros(sz, sx);
                im_names = [""];
                max_ = 0;
                for im = 1:21
                    im_names = [im_names, Ims(im).name];
                    Ims(im).image = Ims(im).image(1:sz, end-sx+1:end);
                    if max(abs(Ims(im).image(:))) > max_
                        max_ = max(abs(Ims(im).image(:))); 
                    end
%                     Ims(im).image = Ims(im).image .* ~global_mask(1:sz, end-sx+1:end, im);
                    E_s = Sens(im).image(1:sz, end-sx+1:end) ./ sigma(1:sz, end-sx+1:end, im);
                    metric = E_s + sqrt(pi/2) ./ (E_s / sqrt(8*pi) + 1).^4;
                    higher = abs(metric) > abs(metric_best);
                    metric_best(higher) = metric(higher);
                    metric_arg(higher) = im;
                end
                metric_arg = reshape(metric_arg, sz, sx);
                Fusion = Ims(19);
                Fusion.image = zeros(sz, sx);
                for im = 1:21
                    Fusion.image(metric_arg == im) = Ims(im).image(metric_arg == im);
                    Ims(im).db_image = 20 * log10(abs(Ims(im).image) ./ max_);
                end

%                 figure(4)
%                 imagesc(Ims(1).x*1e3, flip(-Ims(1).z*1e3), metric_arg)
%                 colormap(jet(20))
%                 c = colorbar;
%                 c.Ticks = 0:21;
%                 c.Limits = [0, 21];
%                 c.TickLabels = im_names;

                for threshold = 1:size(thresholds, 2)
                    % We have made a detection with this threshold.
                    if any(abs(Fusion.image) > thresholds(threshold), 'all')
                        pass_fail(jj+1, threshold) = true;
                    end
                end

%                 max_ = 0;
%                 best_view = 0;
%                 for im = 1:21
%                     temp_im = Sens(im).image;
%                     temp_im(isnan(temp_im)) = 0;
%                     if mean(abs(temp_im)) > max_
%                         max_ = mean(abs(temp_im));
%                         best_view = im;
%                     end
%                 end
%                 best_views = [best_views; Ims(best_view).name];



                %% ROC on Fisher fused data.
%                 load(sprintf("%s big Fisher fusion Test Statistic", matname))
% 
%                 ims_x  = reshape(Fused.x, 1, []);
%                 ims_z  = reshape(Fused.z, 1, []);
%                 [~, sx]  = min(abs(geometry(3).point1(1) - ims_x));
%                 [~, sz]  = min(abs(geometry(4).point1(3) - ims_z));
%                 [~, ex]  = min(abs(geometry(1).point1(1) - ims_x));
%                 [~, ez]  = min(abs(geometry(2).point1(3) - ims_z));
                
%                 subset = Fused.db_image(sz:ez, sx:ex);
%                 exc = sum(subset(:) < thresholds, 1) / ((ex-sx+1)*(ez-sz+1));
%                 plot(thresholds, exc)

%                 for threshold = 1:size(thresholds, 2)
%                     % We have made a detection with this threshold.
%                     if any(Fused.db_image(sz:ez, sx:ex) > thresholds(threshold), 'all')
%                         pass_fail(jj+1, threshold) = true;
%                     end
%                     if and(strcmp(matname, "01"), contains(folders(ii).name, "jig at corner 2"))
%                         figure(4)
%                         points = zeros(ez-sz+1, ex-sx+1);
%                         points(Fused.db_image(sz:ez, sx:ex) < thresholds(threshold)) = 1;
%                         imagesc(points)
%                         colorbar
%                         colormap gray
%                         clim([0, 1])
%                         drawnow
%                         pause(.01)
%                     end

%                     if sum(abs(Fused.db_image(sz:ez, sx:ex)) <= thresholds(threshold), 'all')/((ex-sx+1)*(ez-sz+1)) >= thresholds(threshold)
%                         pass_fail(jj+1, threshold) = true;
%                     end
%                 end

                a = 1;

%% Old from noise quant
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
%                 [~, x_min] = min(abs(Ims(1).x - geometry(3).point1(1)), [], 'all');
%                 [~, x_max] = min(abs(Ims(1).x - geometry(1).point1(1)), [], 'all');
%                 [~, z_min] = min(abs(Ims(1).z - geometry(4).point1(3)), [], 'all');
%                 [~, z_max] = min(abs(Ims(1).z - geometry(2).point1(3)), [], 'all');
%                 x_min = x_min + 1;
%                 x_max = x_max - 1;
%                 z_min = z_min + 1;
%                 z_max = z_max - 1;
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
%             load(sprintf("%s%s", matname, suffix))
%             h = permute(h, [2,1,3]);
%             [X, Z] = meshgrid(acl_x, acl_z);
%             % Assume Rayleigh distribution
%             S = Ims;
%             [~, x_min] = min(abs(S(1).x - acl_x(1)));
%             [~, x_max] = min(abs(S(1).x - acl_x(end)));
%             [~, z_min] = min(abs(S(1).z - acl_z(1)));
%             [~, z_max] = min(abs(S(1).z - acl_z(end)));
%             mask = zeros(size(h));
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
% 
%                 % Update S based on mask.
% %                 S(im).image(squeeze(mask(:, :, im))) = nan;
% %                 S(im).db_image(squeeze(mask(:, :, im))) = nan;
%                 % Fit plane and correct for it
%                 plane = fit([X(:), Z(:)], S(im).db_image(:), 'poly11', 'Exclude', reshape(mask(:, :, im), [], 1));
%                 c = plane.p10 * X + plane.p01 * Z;
%                 S(im).image = abs(S(im).image) ./ 10.^(c / 20);
%                 S(im).db_image = 20 * log10(abs(S(im).image) ./ max(abs(S(im).image(:)), [], 'omitnan'));
%                 % Fit Rayleigh to image
%                 i_c = 0:max(S(im).image(:), [], 'omitnan')/1000:max(S(im).image(:));
% %                 try
%                 rayleigh = fitdist(S(im).image(:), 'Rayleigh');
% %                 catch ME
% %                     break
% %                 end
%                 % Mark everything >99th percentile
% %                 p = icdf(rayleigh, percentile_target);
%                 T = double(S(im).db_image > -20);
%                 mask(:, :, im) = T;
%             end
% 
%             save(sprintf("%s 20db_thresh mask", matname), "mask", "acl_x", "acl_z")

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
%             mask_type = "acf+hics mask";
%             load(sprintf("%s %s", matname, mask_type))
%             for im = 1:21
%                 Masks(im).geometry = geometry;
%                 Masks(idx).matname = matname;
%                 Masks(idx).x = acl_x;
%                 Masks(idx).z = acl_z;
%                 Masks(im).db_image = global_mask(:, :, im);
%             end
% 
%             idx = idx + 1;
            catch ME
                disp(ME.message)
            end

        end

        rate = sum(pass_fail, 1) / size(pass_fail, 1);
        save(sprintf("%s %s pass fail", folders(ii).name, mask_type), "rate", "pass_fail", "thresholds")
%         samples = samples - 100;
    end
end

% mask_type = "";
%% Use pass/fail to plot ROC.
% We have 100 thresholds.
nd = zeros(size(thresholds));
de = zeros(size(thresholds));
for ii = 1:size(folders, 1)
    if any(strcmp(folders(ii).name, stat_folders))
        thisdir = fullfile(dirname, folders(ii).name, 'TFMs', 'Relative coords');
        cd(thisdir)
        load(sprintf("%s %s pass fail", folders(ii).name, mask_type))
        if contains(folders(ii).name, "non-defective")
            nd = nd + sum(pass_fail, 1);
        else 
            de = de + sum(pass_fail, 1);
        end
    end
end
de = de / 220;
nd = nd / 220;
figure(1)
plot(thresholds, nd)
xlabel("Threshold")
ylabel("PFA")
title("Pristine")
figure(2)
plot(thresholds, de)
xlabel("Threshold")
ylabel("POD")
title("Defective")
figure(3)
plot(nd, de)
xlabel("PFA")
ylabel("POD")
title("best view")
save(sprintf("%s roc", mask_type)

Mask = Ims;
for im = 1:21
    Mask(im).db_image = global_mask(:, :, im);
end
Mask = rmfield(Mask, "plotExtras");
fn_image_from_mat(Mask)
grp = get(get(gcf, 'Children'), 'Children');
for im = 2:22
    grp(im).CLim = [0, 1];
    colormap gray
end

% names = [];
% for im = 1:21
%     names = [names; Ims(im).name];
% end
% 
% out = categorical(best_views, names);
% figure(4)
% histogram(out, 'Normalization', 'pdf')

% Want to redefine acl_x and acl_z so that they are all on the same base.
% Subtract S2.point1.

mask_type = "acf mask";
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
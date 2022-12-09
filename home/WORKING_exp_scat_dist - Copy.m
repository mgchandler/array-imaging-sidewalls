close all
clear
clc

dirname = "C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Experimental Data\2022\Amplitude Distribution Study";
% dirname = "/user/work/mc16535/complex_distribution/data";
% addpath(genpath("/user/work/mc16535/complex_distribution"))
cd(dirname)

yaml_name = "exp_scat_dist.yml";
% yaml_name = "artefact_only.yml";
yaml_options = yaml.loadFile(yaml_name);

b_or_s = "sml";

N = 1000;

% % Analytical
% v_L = 6317.01;
% v_S = 3110.28;
% FE
v_L = 6334.87;
v_S = 3113.59;
% % Exp
% v_L = 2913;
% v_S = 3229;

for kk = 1:size(yaml_options.mesh.geom.z, 2)
    yaml_options.mesh.geom.z{kk} = -yaml_options.mesh.geom.z{kk};
end
for kk = 1:size(yaml_options.mesh.scat.z, 2)
    yaml_options.mesh.scat.z{kk} = -yaml_options.mesh.scat.z{kk};
end
yaml_options.mesh.scat = fn_scat_info( ...
    yaml_options.mesh.scat.type, ...
    yaml_options.mesh.scat.x{1}, ...
    0, ...
    yaml_options.mesh.scat.z{1}, ...
    yaml_options.mesh.scat.r{1}, ...
    v_L/yaml_options.probe.freq, ...
    v_S/yaml_options.probe.freq, ...
    deg2rad(0), ...
    'ang_pts_over_2pi', 120 ...
);
yaml_options.model.time_it = false;
% yaml_options.model.image_range = [yaml_options.mesh.scat.x - 3e-3, yaml_options.mesh.scat.x + 3e-3, yaml_options.mesh.scat.z - 3e-3, yaml_options.mesh.scat.z + 3e-3];
yaml_options.model.pixel = .2e-3;
yaml_options.mesh.geom.n_pts = 2000;
yaml_options.mesh.n_per_wl = 0;
yaml_options.model.savepath = "";%"C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\array-imaging-sidewalls\home-output\scat vs artefact distribution\RT Simulation\Random";

%% Run through data sets
folders = dir(dirname);
ff = 1;
dd = 1;
if strcmp(b_or_s, "big")
    data = zeros(4, 100, 21, 47*68);
else
    data = zeros(4, 100, 21, 31^2);
end
stat_folders = ["defective - jig at corner - 0mm offset", "defective - jig at corner - 1mm offset", "non-defective - jig at corner - 0mm offset", "non-defective - jig at corner - 1mm offset"];
for ii = 1:length(folders)
    if and(folders(ii).isdir, and(~strcmp(folders(ii).name, '.'), ~strcmp(folders(ii).name, '..')))
        thisdir = fullfile(dirname, folders(ii).name, 'TFMs', 'Global coords');
        yaml_options.model.savepath = fullfile(thisdir, 'TFMs');
        for jj = 0:99
            try
    %             cd(thisdir)
%                 matname = sprintf("%02d", jj);
%                 load(fullfile(thisdir, matname))

%                 freq = [0:length(exp_data.time)-1] / exp_data.time(end);
%                 yaml_options.data.time = exp_data.time - 5e-7;
%                 yaml_options.data.data = ifft(2 * fn_hanning(length(exp_data.time), yaml_options.probe.freq/max(freq), yaml_options.probe.freq/max(freq)) .* fft(exp_data.time_data));
%                 yaml_options.model.savename = sprintf("%s %s", matname, b_or_s);
%                 geom = fn_size_geometry_from_fmc(exp_data, yaml_options.mesh.geom.z{4}, yaml_options.mesh.geom.x{1}, true, fullfile(thisdir, 'TFMs', sprintf('%s sizing.mat', yaml_options.model.savename)));
%                 if strcmp(b_or_s, "sml")
%                     yaml_options.model.image_range = [geom(1).point1(1) - 18e-3, geom(1).point1(1) - 12e-3, geom(4).point1(3) + 7e-3, geom(4).point1(3) + 13e-3];
%                 end
%                 model_options = fn_default_model_options(yaml_options);
%                 model_options.mesh.geom.geometry = geom;
%                 disp(model_options.model.savepath)
%                 disp(model_options.model.savename)
%                 [Ims, ~, ~] = fn_tfm(model_options);


                if any(strcmp(folders(ii).name, stat_folders))
                load(fullfile(thisdir, sprintf("%02d sml.mat", jj)))
                    for im = 1:21
    %                 data(ff, jj+1, im, :) = reshape(Ims(im).image(37-4:37+4, 53-4:53+4), [], 1);
                        data(dd, jj+1, im, :) = reshape(Ims(im).image(:), [], 1);
                    end
                end
            catch ME
                continue
            end
        end
        % Histogram it
%         hist(abs(squeeze(data(:, 12, 480))))
        ff = ff + 1;
        if any(strcmp(folders(ii).name, stat_folders))
            dd = dd + 1;
        end
    end
end

de = zeros(2*size(data, 2), size(data, 3), size(data, 4));
nd = zeros(2*size(data, 2), size(data, 3), size(data, 4));
de(1:size(data, 2), :, :)                 = data(1, :, :, :);
de(size(data, 2)+1:2*size(data, 2), :, :) = data(2, :, :, :);
nd(1:size(data, 2), :, :)                 = data(3, :, :, :);
nd(size(data, 2)+1:2*size(data, 2), :, :) = data(4, :, :, :);

N_bins = 15;
edges = linspace(min(abs([de(:); nd(:)])), max(abs([de(:); nd(:)])), N_bins+1);
% Need to loop as histcounts treats the ndarray A as a vector A(:)
de_counts = zeros(N_bins, size(data, 3), size(data, 4));
nd_counts = zeros(N_bins, size(data, 3), size(data, 4));
nd_rician = zeros(2, size(data, 3), size(data, 4));
for pt = 1:size(de, 3)
    for im = 1:size(de, 2)
        de_counts(:, im, pt) = histcounts(abs(de(:, im, pt)), edges, 'Normalization', 'pdf');
        nd_counts(:, im, pt) = histcounts(abs(nd(:, im, pt)), edges, 'Normalization', 'pdf');
%         if im==21
%             a=1;
%         end
        try
            fitted = fitdist(abs(nd(:, im, pt)), 'Rician');
        catch ME
            clear fitted
            fitted.s = nan;
            fitted.sigma = nan;
        end
        nd_rician(:, im, pt) = [fitted.s, fitted.sigma];
        
%         figure
%         bar((edges(1:end-1)+edges(2:end))/2, nd_counts(:, im, pt))
%         hold on
%         y = linspace(min(abs([de(:); nd(:)])), max(abs([de(:); nd(:)])), 1000);
%         plot(y, rice(y, fitted.s, fitted.sigma))
    end
end

% Sigma = Ims(1);
% Nu = Ims(1);
% % for im = 1:21
% %     Nu(im).db_image = reshape(nd_rician(1, im, :), 31, 31);
% %     Sigma(im).db_image = reshape(nd_rician(2, im, :), 31, 31);
% % end
% Nu.db_image = reshape(nd_rician(1, 1, :), 31, 31);
% Sigma.db_image = reshape(nd_rician(2, 1, :), 31, 31);
% fn_image_from_mat(Sigma)
% grp = get(get(gcf, 'Children'), 'Children');
% grp(2).CLim = [0, max(Nu(23-im).db_image(:))];%[0, 12.7];
% grp(2).XLim = [Nu.x(1)*1e3, Nu(1).x(end)*1e3];
% grp(2).YLim = [Nu.z(1)*1e3, Nu(1).z(end)*1e3];
% grp(1).Label.String = "\nu";
% grp(1).Label.String = "\sigma";
% grp = get(get(gcf, 'Children'), 'Children');
% for im = 2:22
%     grp(im).CLim = [0, max(Nu(23-im).db_image(:))];%[0, 12.7];
%     grp(im).XLim = [Nu(im-1).x(1)*1e3, Nu(im-1).x(end)*1e3];
%     grp(im).YLim = [Nu(im-1).z(1)*1e3, Nu(im-1).z(end)*1e3];
% end
% grp(1).Label.String = "\nu";
% fn_image_from_mat(Sigma)
% grp = get(get(gcf, 'Children'), 'Children');
% for im = 2:22
%     grp(im).CLim = [0, max(Sigma(23-im).db_image(:))];%[0, 3.2];
%     grp(im).XLim = [Sigma(im-1).x(1)*1e3, Sigma(im-1).x(end)*1e3];
%     grp(im).YLim = [Sigma(im-1).z(1)*1e3, Sigma(im-1).z(end)*1e3];
% end
% grp(1).Label.String = "\sigma";

%% Fisher fusion method
p_threshold = 1;
for ii = 1:length(folders)
%     if and(folders(ii).isdir, and(~strcmp(folders(ii).name, '.'), and(~strcmp(folders(ii).name, '..'), any(strcmp(folders(ii).name, stat_folders)))))
    if and(folders(ii).isdir, and(~strcmp(folders(ii).name, '.'), ~strcmp(folders(ii).name, '..')))
        thisdir = fullfile(dirname, folders(ii).name);
        yaml_options.model.savepath = fullpath(thisdir, 'TFMs');
        for jj = 0:99
            try
%                 cd(thisdir)
                matname = sprintf("%02d %s", jj, b_or_s);
                load(fullfile(thisdir, matname))
                data = zeros(size(Ims, 1), size(Ims(1).image, 1), size(Ims(1).image, 2));
                for im = 1:21
                    data(im, :, :) = abs(Ims(im).image);
                end
                probs = zeros(size(data));
                for im = 1:21
                    pt = 1;
                    for zpt = 1:size(data, 3)
                        for xpt = 1:size(data, 2)
                            probs(im, xpt, zpt) = 1 - cdf('Rician', data(im, xpt, zpt), nd_rician(1, im, pt), nd_rician(2, im, pt));
                            pt = pt + 1;
                        end
                    end
                end
                Probs = repmat(struct('image', zeros(size(probs, 2), size(probs, 3))), size(probs, 1), 1);
                for im = 1:21
                    Probs(im).db_image = squeeze(probs(im, :, :));
                    Probs(im).x = Ims(1).x;
                    Probs(im).z = Ims(1).z;
                    Probs(im).name = Ims(im).name;
                end
                fn_image_from_mat(Probs)
                grp = get(get(gcf, 'Children'), 'Children');
                for im = 2:22
                    grp(im).CLim = [0, p_threshold];
                end
                grp(1).Label.String = "p";
                savefig(fullfile(thisdir, 'TFMs', sprintf("%s View-wise Probability.fig", matname)))
                yaml_options.model.savename = sprintf("%s Fisher fusion", matname);
                fn_fisher_fusion(yaml_options, Probs);
                close all
            catch ME
                continue
            end 
        end
    end
end

%% Likelihood ratio method
for ii = 1:length(folders)
%     if and(folders(ii).isdir, and(~strcmp(folders(ii).name, '.'), and(~strcmp(folders(ii).name, '..'), any(strcmp(folders(ii).name, stat_folders)))))
    if and(folders(ii).isdir, and(~strcmp(folders(ii).name, '.'), ~strcmp(folders(ii).name, '..')))
        thisdir = fullfile(dirname, folders(ii).name, "TFMs");
        yaml_options.model.savepath = thisdir;
        for jj = 0:99
            try
%                 cd(thisdir)
                matname = sprintf("%02d %s", jj, b_or_s);
                load(fullfile(thisdir, sprintf("%s.mat", matname)))
                data = zeros(size(Ims, 1), size(Ims(1).image, 1), size(Ims(1).image, 2));
                for im = 1:21
                    data(im, :, :) = abs(Ims(im).image);
                end
                T_stat = repmat(struct('image', zeros(size(data, 2), size(data, 3))), 1, 1);
                Likelihood = repmat(struct('image', zeros(size(data, 2), size(data, 3))), size(Ims, 1), 1);
                T_stat.x = Ims(1).x;
                T_stat.z = Ims(1).z;
                T_stat.name = "Linear Signal Likelihood";
                likelihood = zeros(size(Ims, 1), size(data, 2), size(data, 3));
                for im = 1:21
                    pt = 1;
                    for zpt = 1:size(data, 3)
                        for xpt = 1:size(data, 2)
                            xi = abs(data(im, xpt, zpt));
                            s = nd_rician(1, im, pt);
                            sigma = nd_rician(2, im, pt);
                            likelihood(im, xpt, zpt) = (xi/sigma)^2 - 2*log(besseli(0, xi*s/sigma^2));
                            T_stat.image(xpt, zpt) = T_stat.image(xpt, zpt) + (xi/sigma)^2 - 2*log(besseli(0, xi*s/sigma^2));
                            pt = pt + 1;
                        end
                    end
                end
                for im = 1:size(Ims, 1)
                    Likelihood(im).db_image = squeeze(likelihood(im, :, :));
                    Likelihood(im).x = Ims(1).x;
                    Likelihood(im).z = Ims(1).z;
                    Likelihood(im).name = Ims(im).name;
                end
                fn_image_from_mat(Likelihood)
                grp = get(get(gcf, 'Children'), 'Children');
                max_ = 0;
                for im = 1:21
                    if max(Likelihood(im).db_image(:)) > max_
                        max_ = max(Likelihood(im).db_image(:));
                    end
                end
                for im = 2:22
                    grp(im).CLim = [0, max_];
                end
                grp(1).Label.String = "Likelihood";
                savefig(fullfile(thisdir, sprintf("%s View-wise Likelihood.fig", matname)))
                
                T_stat.db_image = T_stat.image;
                fn_image_from_mat(T_stat)
                grp = get(get(gcf, 'Children'), 'Children');
                grp(2).CLim = [nanmin(T_stat.db_image(:)), nanmax(T_stat.db_image(:))];
                grp(2).Label.String = "Likelihood";
                savefig(fullfile(thisdir, sprintf("%s Linear Likelihood.fig", matname)))
                close all
            catch ME
                continue
            end 
        end
    end
end



%% Support fns

function y = rice(x, s, sigma)
y = besseli(0, x*s/sigma^2) .* x ./ sigma^2 .* exp(-((x.^2 + s^2) ./ (2*sigma^2)));
end
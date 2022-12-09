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
load_ims = true;
load_in_data = false;
plot_everything = true;

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
    xsize = 70e-3;
    zsize = 50e-3;
    pixel = .5e-3;
else
    xsize = 3e-3;
    zsize = 3e-3;
    pixel = .1e-3;
end
data = zeros(4, 100, 21, (xsize/pixel+1)*(zsize/pixel+1));

stat_folders = ["defective - jig at corner - 0mm offset", "defective - jig at corner - 1mm offset", "non-defective - jig at corner - 0mm offset", "non-defective - jig at corner - 1mm offset"];
for ii = 1:length(folders)
    if and(folders(ii).isdir, and(~strcmp(folders(ii).name, '.'), ~strcmp(folders(ii).name, '..')))
        thisdir = fullfile(dirname, folders(ii).name, 'TFMs', 'Relative coords');
        yaml_options.model.savepath = thisdir;
        for jj = 0:99
            try
                if ~load_ims
                %% Do the imaging
                    cd(thisdir)
                    matname = sprintf("%02d", jj);
                    load(fullfile(thisdir, matname))
                    freq = [0:length(exp_data.time)-1] / exp_data.time(end);
                    yaml_options.data.time = exp_data.time - 5e-7;
                    yaml_options.data.data = ifft(2 * fn_hanning(length(exp_data.time), yaml_options.probe.freq/max(freq), yaml_options.probe.freq/max(freq)) .* fft(exp_data.time_data));
                    yaml_options.model.savename = sprintf("%s %s", matname, b_or_s);
                    geom = fn_size_geometry_from_fmc(exp_data, yaml_options.mesh.geom.z{4}, yaml_options.mesh.geom.x{1}, true, fullfile(thisdir, 'TFMs', sprintf('%s sizing.mat', yaml_options.model.savename)));
                    if strcmp(b_or_s, "sml")
                        yaml_options.model.image_range = [geom(1).point1(1) - 18e-3, geom(1).point1(1) - 12e-3, geom(4).point1(3) + 7e-3, geom(4).point1(3) + 13e-3];
                    end
                    model_options = fn_default_model_options(yaml_options);
                    model_options.mesh.geom.geometry = geom;
                    disp(model_options.model.savepath)
                    disp(model_options.model.savename)
                    [Ims, ~, ~] = fn_tfm(model_options);
                end

                %% Store the data
                if any(strcmp(folders(ii).name, stat_folders))
                    if load_ims
                        load(fullfile(thisdir, sprintf("%02d %s.mat", jj, b_or_s)))
                    end
                    for im = 1:21
    %                 data(ff, jj+1, im, :) = reshape(Ims(im).image(37-4:37+4, 53-4:53+4), [], 1);
                        data(dd, jj+1, im, :) = reshape(Ims(im).image(:), [], 1);
                    end
                end
            catch ME
                continue
            end
        end
        ff = ff + 1;
        if any(strcmp(folders(ii).name, stat_folders))
            dd = dd + 1;
        end
    end
end

% Unpack data
de = zeros(2*size(data, 2), size(data, 3), size(data, 4));
nd = zeros(2*size(data, 2), size(data, 3), size(data, 4));
de(1:size(data, 2), :, :)                 = data(1, :, :, :);
de(size(data, 2)+1:2*size(data, 2), :, :) = data(2, :, :, :);
nd(1:size(data, 2), :, :)                 = data(3, :, :, :);
nd(size(data, 2)+1:2*size(data, 2), :, :) = data(4, :, :, :);


%% Fit to the distributions
% Initialise distribution storage + goodness of fit storage
ksp = zeros(21, size(Ims(1).image, 1)*size(Ims(1).image, 2));
ksv = zeros(21, size(Ims(1).image, 1)*size(Ims(1).image, 2));
chi2p = zeros(21, size(Ims(1).image, 1)*size(Ims(1).image, 2));
chi2v = zeros(21, size(Ims(1).image, 1)*size(Ims(1).image, 2));
nd_rician = zeros(2, size(data, 3), size(data, 4));
if load_in_data
    load(sprintf("%s rician params.mat", b_or_s))
else
    for pt = 1:size(de, 3)
        for im = 1:size(de, 2)
            try
                %% Do the fitting
%                 fitted = fitdist(abs(nd(:, im, pt)), 'Rician');
                
                nu_temp = mean(nd(:, im, pt));
                temp_nd = real(nd(:, im, pt).*exp(-1j*angle(nu_temp)));
                sigma   = std(temp_nd) / sqrt(2);
                nu = mean(abs(nd(:, im, pt)));
                
%                 figure
%                 hold on; axis equal
%                 scatter(real(nd(:, im, pt)), imag(nd(:, im, pt)), '.', 'MarkerEdgeColor', "#0072BD", 'MarkerFaceColor', "#0072BD")
%                 scatter(real(nd(:, im, pt).*exp(-1j*angle(nu_temp))), imag(nd(:, im, pt).*exp(-1j*angle(nu_temp))), '.', 'MarkerEdgeColor', "#77AC30", 'MarkerFaceColor', "#0072BD")
%                 scatter(real(nu_temp), imag(nu_temp), 'o', 'filled', 'MarkerEdgeColor', "#D95319", 'MarkerFaceColor', "#D95319")
%                 scatter(nu, 0, 'ro', 'filled')
%                 plot([0, 0], get(gca, 'YLim'), 'Color', [.5, .5, .5], 'LineStyle', '--')
%                 plot(get(gca, 'XLim'), [0, 0], 'Color', [.5, .5, .5], 'LineStyle', '--')
%                 xlabel("Real(x)")
%                 ylabel("Imag(x)")
%                 figure
%                 hold on
%                 histogram(abs(nd(:, im, pt)), 'Normalization', 'pdf')
%                 x = [0:0.01:12];
%                 plot(x, rice(x, fitted.s, fitted.sigma), 'LineWidth', 2)
%                 plot(x, rice(x, abs(mean(nd(:, im, pt))), std(nd(:, im, pt))/sqrt(2)), 'LineWidth', 2)
%                 plot(x, rice(x, nu, sigma), 'LineWidth', 2)
%                 legend('LL-LL Middle of region', 'fitdist', 'mean, std', 'e^{-i\theta} mean, std')
                
                % Evaluate the goodness of fit
                [~, p, s] = kstest(abs(nd(:, im, pt)), 'CDF', makedist("Rician", "s", nu, "sigma", sigma));
                ksp(im, pt) = p;
                ksv(im, pt) = s;
                [~, p, s] = chi2gof(abs(nd(:, im, pt)), 'CDF', makedist("Rician", "s", nu, "sigma", sigma));
                chi2p(im, pt) = p;
                chi2v(im, pt) = s.chi2stat;
            catch ME
                %% If fit can't be found, store as nan
                disp(ME.message)
%                 clear fitted
%                 fitted.s = nan;
%                 fitted.sigma = nan;
                nu = nan;
                sigma = nan;
                ksp(im, pt) = nan;
                ksv(im, pt) = nan;
                chi2p(im, pt) = nan;
                chi2v(im, pt) = nan;
            end
            nd_rician(:, im, pt) = [nu, sigma];
        end
    end
%     save(fullfile(dirname, sprintf("%s rician params.mat", b_or_s)), "nd_rician", "ksp", "ksv", "chi2p", "chi2v")
end

%% Plot fit parameters
if plot_everything
    Sigma = Ims;
    Nu = Ims;
    for im = 1:21
        Nu(im).db_image = reshape(nd_rician(1, im, :), size(Ims(1).image, 1), size(Ims(1).image, 2));
        Sigma(im).db_image = reshape(nd_rician(2, im, :), size(Ims(1).image, 1), size(Ims(1).image, 2));
    end
    fn_image_from_mat(Nu)
    grp = get(get(gcf, 'Children'), 'Children');
    grp(1).Label.String = "\nu";
    grp = get(get(gcf, 'Children'), 'Children');
    for im = 2:22
        grp(im).CLim = [0, 15];%[0, max(Nu(23-im).db_image(:))];%
        grp(im).XLim = [Nu(im-1).x(1)*1e3, Nu(im-1).x(end)*1e3];
        grp(im).YLim = [Nu(im-1).z(1)*1e3, Nu(im-1).z(end)*1e3];
    end
    grp(1).Label.String = "\nu";
    savefig(fullfile(dirname, sprintf("Nu %s Non-Defective Relative Coords estimated.fig", b_or_s)))
    fn_image_from_mat(Sigma)
    grp = get(get(gcf, 'Children'), 'Children');
    for im = 2:22
        grp(im).CLim = [0, 3];%[0, max(Sigma(23-im).db_image(:))];%
        grp(im).XLim = [Sigma(im-1).x(1)*1e3, Sigma(im-1).x(end)*1e3];
        grp(im).YLim = [Sigma(im-1).z(1)*1e3, Sigma(im-1).z(end)*1e3];
    end
    grp(1).Label.String = "\sigma";
    savefig(fullfile(dirname, sprintf("Sigma %s Non-Defective Relative Coords estimated.fig", b_or_s)))

%% Plot goodness of fit parameters
    KS = Ims;
    Chi = Ims;
    KSv = Ims;
    Chiv = Ims;
    for im = 1:21
        KS(im).db_image = reshape(ksp(im, :), size(Ims(1).image, 1), size(Ims(1).image, 2));
        Chi(im).db_image = reshape(chi2p(im, :), size(Ims(1).image, 1), size(Ims(1).image, 2));
        KSv(im).db_image = reshape(ksv(im, :), size(Ims(1).image, 1), size(Ims(1).image, 2));
        Chiv(im).db_image = reshape(chi2v(im, :), size(Ims(1).image, 1), size(Ims(1).image, 2));
    end
    fn_image_from_mat(KS)
    grp = get(get(gcf, 'Children'), 'Children');
    for im = 2:22
        grp(im).CLim = [0, .05];%[0, max(Sigma(23-im).db_image(:))];%
        grp(im).XLim = [KS(im-1).x(1)*1e3, KS(im-1).x(end)*1e3];
        grp(im).YLim = [KS(im-1).z(1)*1e3, KS(im-1).z(end)*1e3];
    end
    grp(1).Label.String = "p-val from KS Statistic";
    savefig(fullfile(dirname, sprintf("KS p-val %s estimated.fig", b_or_s)))
    fn_image_from_mat(Chi)
    grp = get(get(gcf, 'Children'), 'Children');
    for im = 2:22
        grp(im).CLim = [0, .05];%[0, max(Sigma(23-im).db_image(:))];%
        grp(im).XLim = [Chi(im-1).x(1)*1e3, Chi(im-1).x(end)*1e3];
        grp(im).YLim = [Chi(im-1).z(1)*1e3, Chi(im-1).z(end)*1e3];
    end
    grp(1).Label.String = "p-val from \chi^2 test";
    savefig(fullfile(dirname, sprintf("Chi2 p-val %s estimated.fig", b_or_s)))
    fn_image_from_mat(KSv)
    grp = get(get(gcf, 'Children'), 'Children');
    for im = 2:22
        grp(im).CLim = [0, .3];%[0, max(Sigma(23-im).db_image(:))];%
        grp(im).XLim = [KSv(im-1).x(1)*1e3, KSv(im-1).x(end)*1e3];
        grp(im).YLim = [KSv(im-1).z(1)*1e3, KSv(im-1).z(end)*1e3];
    end
    grp(1).Label.String = "KS Statistic";
    savefig(fullfile(dirname, sprintf("KS Statistic %s estimated.fig", b_or_s)))
    fn_image_from_mat(Chiv)
    grp = get(get(gcf, 'Children'), 'Children');
    for im = 2:22
        grp(im).CLim = [0, 85];%[0, max(Sigma(23-im).db_image(:))];%
        grp(im).XLim = [Chiv(im-1).x(1)*1e3, Chiv(im-1).x(end)*1e3];
        grp(im).YLim = [Chiv(im-1).z(1)*1e3, Chiv(im-1).z(end)*1e3];
    end
    grp(1).Label.String = "\chi^2";
    savefig(fullfile(dirname, sprintf("Chi2 statistic %s estimated.fig", b_or_s)))
end



%% Fisher fusion method
p_threshold = 1;
for ii = 1:length(folders)
%     if and(folders(ii).isdir, and(~strcmp(folders(ii).name, '.'), and(~strcmp(folders(ii).name, '..'), any(strcmp(folders(ii).name, stat_folders)))))
    if and(folders(ii).isdir, and(~strcmp(folders(ii).name, '.'), ~strcmp(folders(ii).name, '..')))
        thisdir = fullfile(dirname, folders(ii).name, 'TFMs', 'Relative coords');
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
                probs = zeros(size(data));
                chinum = zeros(size(data));
                for im = 1:21
                    pt = 1;
                    for zpt = 1:size(data, 3)
                        for xpt = 1:size(data, 2)
                            if ksp(im, pt) > 0.05
                                probs(im, xpt, zpt) = 1 - cdf('Rician', data(im, xpt, zpt), nd_rician(1, im, pt), nd_rician(2, im, pt));
                                chinum(im, xpt, zpt) = 1;
                            end
                            pt = pt + 1;
                        end
                    end
                end
                Probs = repmat(struct('image', ones(size(probs, 2), size(probs, 3))), size(probs, 1), 1);
                for im = 1:21
                    Probs(im).db_image = squeeze(probs(im, :, :));
                    Probs(im).image = squeeze(chinum(im, :, :));
                    Probs(im).x = Ims(1).x;
                    Probs(im).z = Ims(1).z;
                    Probs(im).name = Ims(im).name;
                    Probs(im).plotExtras = Ims(1).plotExtras;
                end
                fn_image_from_mat(Probs)
                grp = get(get(gcf, 'Children'), 'Children');
                for im = 2:22
                    grp(im).CLim = [0, p_threshold];
                end
                grp(1).Label.String = "p";
                savefig(fullfile(thisdir, sprintf("%s View-wise Probability.fig", matname)))
                yaml_options.model.savename = sprintf("%s Fisher fusion", matname);
                fn_fisher_fusion(yaml_options, Probs);
                close all
            catch ME
                continue
            end 
        end
    end
end

%% Linear Likelihood ratio method
for ii = 1:length(folders)
%     if and(folders(ii).isdir, and(~strcmp(folders(ii).name, '.'), and(~strcmp(folders(ii).name, '..'), any(strcmp(folders(ii).name, stat_folders)))))
    if and(folders(ii).isdir, and(~strcmp(folders(ii).name, '.'), ~strcmp(folders(ii).name, '..')))
        thisdir = fullfile(dirname, folders(ii).name, "TFMs", 'Relative coords');
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
                Lin_Likelihood = repmat(struct('image', zeros(size(data, 2), size(data, 3))), size(Ims, 1), 1);
                T_stat.x = Ims(1).x;
                T_stat.z = Ims(1).z;
                T_stat.name = "Linear Signal Likelihood";
                T_stat.plotExras = Ims(1).plotExtras;
                lin_likelihood = zeros(size(Ims, 1), size(data, 2), size(data, 3));
                for im = 1:21
                    pt = 1;
                    for zpt = 1:size(data, 3)
                        for xpt = 1:size(data, 2)
                            xi = abs(data(im, xpt, zpt));
                            s = nd_rician(1, im, pt);
                            sigma = nd_rician(2, im, pt);
                            lin_likelihood(im, xpt, zpt) = (xi/sigma)^2 - 2*(sqrt(4 + (xi*s/sigma^2)^2) - 2);%log(besseli(0, xi*s/sigma^2));
                            T_stat.image(xpt, zpt) = T_stat.image(xpt, zpt) + (xi/sigma)^2 - 2*log(besseli(0, xi*s/sigma^2));
                            pt = pt + 1;
                        end
                    end
                end
                for im = 1:size(Ims, 1)
                    Lin_Likelihood(im).db_image = squeeze(lin_likelihood(im, :, :));
                    Lin_Likelihood(im).x = Ims(1).x;
                    Lin_Likelihood(im).z = Ims(1).z;
                    Lin_Likelihood(im).name = Ims(im).name;
                    Lin_Likelihood(im).plotExtras = Im(im).plotExtras;
                end
                fn_image_from_mat(Lin_Likelihood)
                grp = get(get(gcf, 'Children'), 'Children');
                max_ = 0;
                for im = 1:21
                    if max(Lin_Likelihood(im).db_image(:)) > max_
                        max_ = max(Lin_Likelihood(im).db_image(:));
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
                grp(1).Label.String = "Likelihood";
                savefig(fullfile(thisdir, sprintf("%s Linear Likelihood.fig", matname)))
                saveas(gcf, fullfile(thisdir, sprintf("%s Linear Likelihood.png", matname)))
                close all
            catch ME
                disp(ME.message)
                continue
            end 
        end
    end
end



%% Support fns
function y = rice(x, s, sigma)
y = besseli(0, x*s/sigma^2) .* x ./ sigma^2 .* exp(-((x.^2 + s^2) ./ (2*sigma^2)));
end
clear all
close all
clc

dirname = "C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Experimental Data\2022\L-shape FMCs with 1mm SDHs\22.09.12 - test + ind sample\FE";
cd(dirname)

yaml_name = "a_0mm_d_1.yml";
% yaml_name = "a 0mm - d 1.yml";
fmc_name  = "a_0mm_d_1_1_BP.mat";
% fmc_name  = "a 0mm - d 1.mat";
bln_name  = "a_0mm_d_1_b_BP.mat";
% bln_name  = "a 0mm - c 1.mat";

yaml_options = yaml.loadFile(yaml_name);
load(bln_name);
threshold = .5e-12;
half_dist = 20;
% fmc_mask.data = ones(size(data));
% fmc_mask.data(abs(data) > threshold) = 0;
% new_mask = ones(size(data));
% for ii = 1:4096
% %     fmc_mask.data(:, ii) = smooth(fmc_mask.data(:, ii), 20);
%     for jj = 1:size(data, 1)
%         if ~fmc_mask.data(jj, ii)
%             new_mask(max(1, jj-half_dist):min(size(data, 1), jj+half_dist), ii) = 0;
%         end
%     end
% end
% fmc_mask.data = new_mask;
% fmc_mask.time = time(:, 1);
bln_data = data;
bln_time = time;
load(fmc_name);
% data = data - bln_data;

% % Analytical
% v_L = 6317.01;
% v_S = 3110.28;
% FE
v_L = 6334.87;
v_S = 3113.59;
% % Exp
% v_L = 2913;
% v_S = 3229;

%% Additional model_options
for kk = 1:size(yaml_options.mesh.geom.z, 2)
    yaml_options.mesh.geom.z{kk} = -yaml_options.mesh.geom.z{kk};
end
for kk = 1:size(yaml_options.mesh.scat.z, 2)
    yaml_options.mesh.scat.z{kk} = -yaml_options.mesh.scat.z{kk};
end
yaml_options.mesh.n_per_wl = 0;
% yaml_options.model.image_range = [30e-3, 48e-3, 26e-3, 45e-3];
% yaml_options.model.image_range = [yaml_options.mesh.scat.x{1}-2.5e-3, yaml_options.mesh.scat.x{1}+2.5e-3, yaml_options.mesh.scat.z{1}-2.5e-3, yaml_options.mesh.scat.z{1}+2.5e-3];
yaml_options.model.pixel = 1.0e-3;
yaml_options.model.wall_for_imaging = "S1";
yaml_options.model.max_no_reflections = 1;
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
...    'fmc_mask', fmc_mask ...
);
% yaml_options.mesh.scat = fn_scat_info( ...
%     "image", ...
%     yaml_options.mesh.scat.x{1}, ...
%     0, ...
%     yaml_options.mesh.scat.z{1} ...
% );

yaml_options.model.savepath = "C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Experimental Data\2022\L-shape FMCs with 1mm SDHs\22.09.12 - test + ind sample\FE";
yaml_options.model.savename = "a_0mm_d_1_1_sens_unmasked";

% freq = [0:length(exp_data.time)-1] / exp_data.time(end);
% time_data = ifft(2 * fn_hanning(length(exp_data.time), yaml_options.probe.freq/max(freq), yaml_options.probe.freq/max(freq)) .* fft(exp_data.time_data));
% 
% yaml_options.data.time = exp_data.time(:, 1) - 5e-7;
% yaml_options.data.data = time_data;

yaml_options.data.time = time(:, 1);
yaml_options.data.data = data;
model_options = fn_default_model_options(yaml_options);

%% Make the RT mask

% % Diffraction
% diffraction_options = model_options;
% diffraction_options.model.image_locs = [diffraction_options.mesh.geom.x{4} + .01e-3, 0.0, diffraction_options.mesh.geom.z{4} - .01e-3];
% diffraction_options.mesh.scat.x = diffraction_options.mesh.geom.x{4} + .01e-3;
% diffraction_options.mesh.scat.z = diffraction_options.mesh.geom.z{4} - .01e-3;
% diffraction_options.model.max_no_reflections = 0;
% diffraction_options.model.model_geom = 0;
% [~, diffraction_Views] = fn_sens(diffraction_options);
% [diff_FMC_time, diff_FMC_data] = fn_simulate_fmc(diffraction_options);

% Geometry
geom_options = model_options;
geom_options.mesh.scat = fn_scat_info("image", geom_options.mesh.scat.x, geom_options.mesh.scat.y, geom_options.mesh.scat.z);
geom_options.model.model_geom = 2;
geom_options.model.max_no_reverberations = 2;
geom_options.model.wall_for_imaging = 'B1';
geom_options.mesh.geom.geometry = fn_make_geometry(0, 0, 5000, [-20e-3, 0, 20e-3], [20e-3, 0, 20e-3], [-20e-3, 0, 40e-3], [20e-3, 0, 40e-3]);
[geom_FMC_time, geom_FMC_data] = fn_simulate_fmc(geom_options);
fn_image_tfm(geom_FMC_time, geom_FMC_data, geom_options);
imagesc(geom_FMC_time, 1:4096, abs(geom_FMC_data.'))

%% Fusion
[Fused_Defect, Weighted_Defect, Ims_Defect] = fn_simple_fusion(model_options);
% [Fused_RT, Weighted_RT, Ims_RT] = fn_fusion_rtmask(model_options);
mask = repmat(struct('fusion_mask', 0), size(Ims_Defect, 1), 1);
for im = 1:size(Ims_Defect, 1)
    mask(im).fusion_mask = -Ims_Defect(im).db_image/model_options.model.db_range;
    mask(im).fusion_mask(mask(im).fusion_mask > 1) = 1;
    mask(im).fusion_mask(mask(im).fusion_mask < .9) = 0;
end
% Fused_Defect.db_image(Fused_Defect.db_image < -40) = -40;
% yaml_options.model.fusion_mask = mask;

fmc_name  = "a 0mm - d 1.mat";
% fmc_name  = "a_0mm_d_1_1_BP.mat";
cd(dirname)
load(fmc_name);

% load(bln_name);
% bln_data = data;
% load(fmc_name);
% data = data - bln_data;

yaml_options.model.savename = "a_0mm_d_1_mask_fused";
yaml_options.data.time = exp_data.time(:, 1) - 5e-7;
freq = [0:length(exp_data.time)-1] / exp_data.time(end);
time_data = ifft(2 * fn_hanning(length(exp_data.time), yaml_options.probe.freq/max(freq), yaml_options.probe.freq/max(freq)) .* fft(exp_data.time_data));
yaml_options.data.data = data;
yaml_options.data.time = time(:, 1);
model_options = fn_default_model_options(yaml_options);
[Fused_Clean, Weighted_Clean, Ims_Clean] = fn_simple_fusion(model_options);
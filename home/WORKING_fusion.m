clear all
close all
clc

dirname = "C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Experimental Data\2022\L-shape FMCs with 1mm SDHs\22.09.12 - test + ind sample\FE";
cd(dirname)

yaml_name = "a_0mm_d_1.yml";
% yaml_name = "a 0mm - d 1.yml";
fmc_name  = "a_0mm_d_1_1_BP.mat";
fmc_name  = "a 0mm - d 1.mat";
bln_name  = "a_0mm_d_1_b_BP.mat";

yaml_options = yaml.loadFile(yaml_name);
load(bln_name);
% bln_data = data;
% load(fmc_name);
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
% yaml_options.model.image_range = [-12.4e-3, 39.9e-3, 0, 45.0e-3];
% yaml_options.model.image_range = [yaml_options.mesh.scat.x{1}-2.5e-3, yaml_options.mesh.scat.x{1}+2.5e-3, yaml_options.mesh.scat.z{1}-2.5e-3, yaml_options.mesh.scat.z{1}+2.5e-3];
yaml_options.model.pixel = 1.0e-3;
yaml_options.model.wall_for_imaging = "S1";
yaml_options.model.number_of_reflections = 1;
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
yaml_options.model.savepath = "C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Experimental Data\2022\L-shape FMCs with 1mm SDHs\22.09.12 - test + ind sample\TFMs";
yaml_options.model.savename = "a_0mm_d_1_b_fused";

% freq = [0:length(exp_data.time)-1] / exp_data.time(end);
% time_data = ifft(2 * fn_hanning(length(exp_data.time), yaml_options.probe.freq/max(freq), yaml_options.probe.freq/max(freq)) .* fft(exp_data.time_data));
% 
% yaml_options.data.time = exp_data.time(:, 1) - 5e-7;
% yaml_options.data.data = time_data;

yaml_options.data.time = time(:, 1);
yaml_options.data.data = data;

%% Fusion
model_options = fn_default_model_options(yaml_options);
[Fused_Defect, Weighted_Defect, Ims_Defect] = fn_simple_fusion(model_options);
mask = repmat(struct('fusion_mask', 0), size(Ims_Defect, 1), 1);
for im = 1:size(Ims_Defect, 1)
    mask(im).fusion_mask = -Ims_Defect(im).db_image/40;
    mask(im).fusion_mask(mask(im).fusion_mask > 1) = 1;
    mask(im).fusion_mask(mask(im).fusion_mask < .8) = 0;
end
% Fused_Defect.db_image(Fused_Defect.db_image < -40) = -40;
yaml_options.model.fusion_mask = mask;

fmc_name  = "a 0mm - d 1.mat";
cd(dirname)
load(fmc_name);
yaml_options.model.savename = "a_0mm_d_1_mask_fused";
yaml_options.data.time = exp_data.time(:, 1) - 5e-7;
freq = [0:length(exp_data.time)-1] / exp_data.time(end);
time_data = ifft(2 * fn_hanning(length(exp_data.time), yaml_options.probe.freq/max(freq), yaml_options.probe.freq/max(freq)) .* fft(exp_data.time_data));
yaml_options.data.data = time_data;
model_options = fn_default_model_options(yaml_options);
[Fused_Clean, Weighted_Clean, Ims_Clean] = fn_simple_fusion(model_options);
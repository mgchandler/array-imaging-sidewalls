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
threshold = 1e-12;
half_dist = 20;
fmc_mask.data = ones(size(data));
fmc_mask.data(abs(data) > threshold) = 0;
new_mask = ones(size(data));
for ii = 1:4096
%     fmc_mask.data(:, ii) = smooth(fmc_mask.data(:, ii), 20);
    for jj = 1:size(data, 1)
        if ~fmc_mask.data(jj, ii)
            new_mask(max(1, jj-half_dist):min(size(data, 1), jj+half_dist), ii) = 0;
        end
    end
end
fmc_mask.data = new_mask;
fmc_mask.time = time(:, 1);
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
    'ang_pts_over_2pi', 120, ...
    'fmc_mask', fmc_mask ...
);
% yaml_options.mesh.scat = fn_scat_info( ...
%     "image", ...
%     yaml_options.mesh.scat.x{1}, ...
%     0, ...
%     yaml_options.mesh.scat.z{1} ...
% );

yaml_options.model.savepath = "C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Experimental Data\2022\L-shape FMCs with 1mm SDHs\22.09.12 - test + ind sample\FE";
yaml_options.model.savename = "a_0mm_d_1_b_maskedsens_fused";

% freq = [0:length(exp_data.time)-1] / exp_data.time(end);
% time_data = ifft(2 * fn_hanning(length(exp_data.time), yaml_options.probe.freq/max(freq), yaml_options.probe.freq/max(freq)) .* fft(exp_data.time_data));
% 
% yaml_options.data.time = exp_data.time(:, 1) - 5e-7;
% yaml_options.data.data = time_data;

yaml_options.data.time = time(:, 1);
yaml_options.data.data = data;

%% Fusion
model_options = fn_default_model_options(yaml_options);

rot_matrix = [cos(model_options.probe.angle) 0 sin(model_options.probe.angle); 0 1 0; -sin(model_options.probe.angle) 0 cos(model_options.probe.angle)];
probe_coords = zeros(3, model_options.probe.num_els);
probe_coords(1, :) = linspace(0, (model_options.probe.num_els - 1) * model_options.probe.separation + model_options.probe.width, model_options.probe.num_els);
probe_coords(1, :) = probe_coords(1, :) - mean(probe_coords(1, :));
probe_coords = probe_coords.' * rot_matrix;
probe_coords(:, 3) = probe_coords(:, 3) - model_options.probe.standoff;
tx_path = fn_path_info( ...
            "L", ...
            "L", ...
            [0], ...
            0, ...
            [model_options.material.v_L], ...
            [model_options.material.couplant_v, model_options.material.v_L, model_options.material.v_S], ...
            [1], ... % index for material identity
            [model_options.material.couplant_density, model_options.material.density], ...
            model_options.probe.freq, ...
            model_options.probe.separation + model_options.probe.width, ...
            model_options.probe.width, ...
            probe_coords, ...
            0 ...
        );
rx_path = fn_path_info( ...
            "L", ...
            "L", ...
            [0], ...
            0, ...
            [model_options.material.v_L], ...
            [model_options.material.couplant_v, model_options.material.v_L, model_options.material.v_S], ...
            [1], ... % index for material identity
            [model_options.material.couplant_density, model_options.material.density], ...
            model_options.probe.freq, ...
            model_options.probe.separation + model_options.probe.width, ...
            model_options.probe.width, ...
            probe_coords, ...
            0 ...
        );
fn_plot_FMC_at_time(data, time(:, 1), tx_path, rx_path, [0, 0, 10e-3; 10e-3, 0, 10e-3], "");

[Fused_Defect, Weighted_Defect, Ims_Defect] = fn_simple_fusion(model_options);
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
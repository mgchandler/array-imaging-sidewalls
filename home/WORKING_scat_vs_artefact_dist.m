close all
clear
clc

% dirname = "C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\array-imaging-sidewalls\home-output\scat vs artefact distribution\RT Simulation\Random";
dirname = "C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\array-imaging-sidewalls\home-output\scat vs artefact distribution\FE Simulation";
cd(dirname)

yaml_name = "scat_vs_artefact_dist.yml";
% yaml_name = "artefact_only.yml";
yaml_options = yaml.loadFile(yaml_name);

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
yaml_options.model.image_range = [yaml_options.mesh.scat.x, yaml_options.mesh.scat.x, yaml_options.mesh.scat.z, yaml_options.mesh.scat.z];
% 
yaml_options.mesh.geom.n_pts = 2000;

yaml_options.model.savepath = "";%"C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\array-imaging-sidewalls\home-output\scat vs artefact distribution\RT Simulation\Random";

%% Isolated geometry and sdh values
cd(sprintf("%s\\scat_only", dirname))
load("scat_only_BP.mat")
yaml_options = yaml.loadFile("scat_only.yml");
yaml_options.data.time = time(:, 1);
yaml_options.data.data = data;
yaml_options.mesh.n_per_wl = 0;
yaml_options.model.image_range = [yaml_options.mesh.scat.x{1}, yaml_options.mesh.scat.x{1}, yaml_options.mesh.scat.z{1}, yaml_options.mesh.scat.z{1}];
model_options = yaml_options;
model_options.model.model_geom = false;
model_options.model.savename = "Perturbation no geom";
model_options = fn_default_model_options(model_options);
[Ims, Views_im, ~] = fn_tfm(model_options);
sdh_only = Ims(1).image;

cd(sprintf("%s\\artefact_only", dirname))
load("artefact_only_BP.mat")
yaml_options = yaml.loadFile("artefact_only.yml");
yaml_options.data.time = time(:, 1);
yaml_options.data.data = data;
yaml_options.mesh.n_per_wl = 0;
yaml_options.model.image_range = [0, 0, yaml_options.mesh.geom.z{2}, yaml_options.mesh.geom.z{2}];
model_options = yaml_options;
model_options.model.model_geom = true;
model_options.mesh.scat.type = "image";
model_options.model.savename = "Perturbation no sdh";
for kk = 1:size(model_options.mesh.geom.z, 2)
    model_options.mesh.geom.z{kk} = -model_options.mesh.geom.z{kk};
end
model_options = fn_default_model_options(model_options);
model_options.mesh.scat.x = 0;
model_options.mesh.scat.z = 19.5e-3;
[Ims, ~, ~] = fn_tfm(model_options);
geom_only = Ims(1).image;

%% Geometry only: perturb back wall location wrt sdh.
% rng('default')
% % pd = makedist('Uniform', 'Lower', -3.5e-3, 'Upper', 0e-3);
% pd = makedist('Normal', 'mu', 2e-3, 'sigma', 1e-3);
% perturb_bw = random(pd, 1, N);
% 
% % model_options = yaml_options;
% % model_options.mesh.scat.x = 0;
% % model_options.mesh.scat.z = 19.5e-3;
% model_options.model.model_geom = true;
% model_options.mesh.scat.type = "image";
% model_options.model.savename = "Perturbation no sdh";
% % for kk = 1:size(model_options.mesh.geom.z, 2)
% %     model_options.mesh.geom.z{kk} = model_options.mesh.scat.z;
% % end
% model_options.model.image_locs = [model_options.mesh.scat.x * ones(size(perturb_bw)); zeros(size(perturb_bw)); model_options.mesh.scat.z+perturb_bw].';
% model_options = fn_default_model_options(model_options);
% perturb_bw = - perturb_bw;
% [Ims, ~, ~] = fn_tfm(model_options);
% bw_only_vals = Ims(1).image;

%% Random perturbation of back wall with SDH

cd(dirname)
yaml_name = "scat_vs_artefact_dist.yml";
% yaml_name = "artefact_only.yml";
yaml_options = yaml.loadFile(yaml_name);
yaml_options.mesh.n_per_wl = 0;
for kk = 1:size(yaml_options.mesh.geom.z, 2)
    yaml_options.mesh.geom.z{kk} = -yaml_options.mesh.geom.z{kk};
end
for kk = 1:size(yaml_options.mesh.scat.z, 2)
    yaml_options.mesh.scat.z{kk} = -yaml_options.mesh.scat.z{kk};
end
yaml_options.model.image_range = [yaml_options.mesh.scat.x{1}, yaml_options.mesh.scat.x{1}, yaml_options.mesh.scat.z{1}, yaml_options.mesh.scat.z{1}];


rng('default')
% pd = makedist('Uniform', 'Lower', -.5e-3, 'Upper', 3e-3);
% pd = makedist('Normal', 'mu', 2e-3, 'sigma', 1e-3);
% perturb = random(pd, 1, N);
perturb = [-0.50, -0.30, -0.20, -0.10, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.50, 3.00];
vals = zeros(size(perturb));

t1 = tic;
for ii = 1:length(perturb)
    p = perturb(ii);
    model_options = yaml_options;
    for kk = 2:size(model_options.mesh.geom.z, 2)-1
        model_options.mesh.geom.z{kk} = model_options.mesh.geom.z{kk} + p;
    end
    cd(sprintf("%s\\perturbation%.2f", dirname, p))
    load(sprintf("perturbation%.2f_BP.mat", p))
    model_options.data.time = time(:, 1);
    model_options.data.data = data;
%     model_options.model.savename = sprintf("Perturbation %.2fmm", p*1e3);
%     disp(model_options.model.savename)
%     t = double(toc(t1));
%     model_options.model.image_range = [min(cell2mat(model_options.mesh.scat.x)), max(cell2mat(model_options.mesh.scat.x)), model_options.scat.geom.z, model_options.scat.geom.z];
    model_options = fn_default_model_options(model_options);
    
    [Ims, ~, ~] = fn_tfm(model_options, Views_im);
    vals(ii) = Ims(1).image;
    
    if mod(ii, round(N/20)) == 0
        step_time = double(toc(t1));
        eta_time = step_time * (N/ii-1);
        if step_time > 3600
            step_str = sprintf("%6.2f h", step_time / 3600);
        elseif step_time > 60
            step_str = sprintf("%6.2f m", step_time / 60);
        else
            step_str = sprintf("%6.2f s", step_time);
        end
        if eta_time > 3600
            eta_str = sprintf("%6.2f h", eta_time / 3600);
        elseif eta_time > 60
            eta_str = sprintf("%6.2f m", eta_time / 60);
        else
            eta_str = sprintf("%6.2f s", eta_time);
        end
        fprintf("itn %6d : p %6.3fmm : Run %s : ETA %s\n", ii, p*1e3, step_str, eta_str)
    end
end

save("random_sampling_fe.mat", "geom_only", "sdh_only", "perturb", "vals")

xlims = [min([real(vals), real(geom_only), real(sdh_only)])-.02, max([real(vals), real(geom_only), real(sdh_only)])+.02];
ylims = [min([imag(vals), imag(geom_only), imag(sdh_only)])-.02, max([imag(vals), imag(geom_only), imag(sdh_only)])+.02];

fig = figure;
h(1) = plot([xlims(1), xlims(2)], [0, 0], 'k','HandleVisibility','off');
hold on, axis equal
h(2) = plot([0, 0], [ylims(1), ylims(2)], 'k','HandleVisibility','off');
h(3) = scatter(real(vals), imag(vals), 'filled', 'CData', (perturb-min(perturb)));
h(4) = scatter(real(sdh_only), imag(sdh_only), 'x', 'MarkerEdgeColor', "#D95319", 'LineWidth', 2);
h(5) = scatter(real(geom_only), imag(geom_only), 'x', 'MarkerEdgeColor', "#77AC30", 'LineWidth', 2);
xlabel("Re\{z\}")
ylabel("Im\{z\}")
xlim(xlims)
ylim(ylims)
legend(h([3,4,5]), "Back wall + SDH", "SDH only", "Back wall only")
colorbar();

% p = get(gca, 'Position');
% superX = [real(vals(1))+0.005,   real(vals(1))];
% superY = [imag(vals(1))-0.005,   imag(vals(1))];
% separX = [real(vals(end)), real(vals(end))];
% separY = [imag(vals(end))-0.005, imag(vals(end))];
% superX = (superX - xlims(1)) ./ (xlims(2) - xlims(1)) + p(1);
% superY = (superY - ylims(1)) ./ (ylims(2) - ylims(1)) + p(2);
% separX = (separX - xlims(1)) ./ (xlims(2) - xlims(1)) + p(1);
% separY = (separY - ylims(1)) ./ (ylims(2) - ylims(1)) + p(2);

% super = annotation('textarrow', superX, superY, 'String', 'Back wall superposed with SDH');
% separ = annotation('textarrow', separX, separY, 'String', 'Back wall separated from SDH');
% super.Parent = fig.CurrentAxes;
% super.String = 'Back wall superposed with SDH';
% separ.Parent = fig.CurrentAxes;
% separ.String = 'Back wall separated from SDH';
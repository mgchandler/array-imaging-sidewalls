close all
clear
clc

dirname = "C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\array-imaging-sidewalls\home-output\scat vs artefact distribution\RT Simulation\Random";
cd(dirname)

yaml_name = "scat_vs_artefact_dist.yml";
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
model_options = yaml_options;
model_options.model.model_geom = false;
model_options.model.savename = "Perturbation no geom";
model_options = fn_default_model_options(model_options);
[Ims, Views_im, ~] = fn_tfm(model_options);
sdh_only = Ims(1).image;

model_options = yaml_options;
model_options.model.model_geom = true;
model_options.mesh.scat.type = "image";
model_options.model.savename = "Perturbation no sdh";
for kk = 1:size(model_options.mesh.geom.z, 2)
    model_options.mesh.geom.z{kk} = model_options.mesh.scat.z;
end
model_options = fn_default_model_options(model_options);
[Ims, ~, ~] = fn_tfm(model_options, Views_im);
geom_only = Ims(1).image;

%% Both sdh and geometry: perturb back wall location wrt sdh.
rng('default')
pd = makedist('Uniform', 'Lower', -.5e-3, 'Upper', 3e-3);
perturb = random(pd, 1, N);
% perturb = -.5e-3:.05e-3:1e-3;
% perturb = [perturb, 1e-3:.1e-3:2e-3];
% perturb = [perturb, 2e-3:.2e-3:3e-3];
perturb = unique(perturb);
vals = zeros(size(perturb));

%% Unique only
model_options = yaml_options;
model_options.model.model_geom = true;
model_options.mesh.scat.type = "image";
model_options.model.savename = "Perturbation no sdh";
for kk = 1:size(model_options.mesh.geom.z, 2)
    model_options.mesh.geom.z{kk} = model_options.mesh.scat.z;
end
model_options.model.image_locs = [model_options.mesh.scat.x * ones(size(perturb)); zeros(size(perturb)); model_options.mesh.scat.z+perturb].';
model_options = fn_default_model_options(model_options);
[Ims, ~, ~] = fn_tfm(model_options);
bw_only_perturb = Ims(1).image;

%%
t1 = tic;
for ii = 1:length(perturb)
    p = perturb(ii);
    model_options = yaml_options;
    for kk = 1:size(model_options.mesh.geom.z, 2)
        model_options.mesh.geom.z{kk} = model_options.mesh.geom.z{kk} + p;
    end
%     model_options.model.savename = sprintf("Perturbation %.2fmm", p*1e3);
%     disp(model_options.model.savename)
%     t = double(toc(t1));
%     model_options.model.image_range = [min(cell2mat(model_options.mesh.geom.x)), max(cell2mat(model_options.mesh.geom.x)), 0, max(cell2mat(model_options.mesh.geom.z))];
    model_options = fn_default_model_options(model_options);
    
    [Ims, ~, ~] = fn_tfm(model_options, Views_im);
    vals(ii) = Ims(1).image;
    
    fprintf("itn %4d : %5.3fmm : ETA %.5fs\n", ii, p*1e3, double(toc(t1))*(N/ii-1))
end

xlims = [min([real(vals), real(geom_only), real(sdh_only)])-.02, max([real(vals), real(geom_only), real(sdh_only)])+.02];
ylims = [min([imag(vals), imag(geom_only), imag(sdh_only)])-.02, max([imag(vals), imag(geom_only), imag(sdh_only)])+.02];

fig = figure;
h(1) = plot([xlims(1), xlims(2)], [0, 0], 'k','HandleVisibility','off');
hold on, axis equal
h(2) = plot([0, 0], [ylims(1), ylims(2)], 'k','HandleVisibility','off');
h(3) = scatter(real(vals), imag(vals), 'filled', 'CData', (perturb-min(perturb))*1e3);
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

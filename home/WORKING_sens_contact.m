clear

v_L = 6320;%6317.0122248907810;
v_S = 3130;%3110.2818131859126;

yaml_options = yaml.loadFile("L_sens.yml");
yaml_options.material.couplant_v = 340.0;
yaml_options.material.couplant_density = 1.2;
for y = 1:size(yaml_options.mesh.geom.y, 2)
    yaml_options.mesh.geom.y{y} = -yaml_options.mesh.geom.y{y};
end
yaml_options.model.boxsize = 1.0e-3;
yaml_options.model.pixel = 1.0e-3;
yaml_options.model.model_geom = 0;
yaml_options.probe.angle = 0.0;
yaml_options.probe.standoff = 0.0;
yaml_options.model.max_no_reflections = 0;
yaml_options.model.savepath = 'C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\array-imaging-sidewalls\array-imaging-sidewalls matlab\output';
yaml_options.model.savename = 'SDH Refl - Valid Paths';
% new_options.geometry = fn_make_geometry(0, 500, ...
% ...%     [new_options.geom_shape.xmax, 0.0, -1e-5], [new_options.geom_shape.xmin, 0.0, -1e-5], ...
%     [new_options.geom_shape.xmax, 0.0, new_options.geom_shape.zmax+1e-5], [new_options.geom_shape.xmin, 0.0, new_options.geom_shape.zmax+1e-5] ...
% );
yaml_options.mesh.sdh.info = fn_scat_info( ...
    "sdh", ...
    .5e-3, ...
    v_L/5e6, ...
    v_S/5e6, ...
    deg2rad(0), ...
    [[32.5e-3, 0.0, 27.5e-3]], ...
    'ang_pts_over_2pi', 120 ...
);
yaml_options.model.wall_for_imaging = "S2";
model_options = fn_default_model_options(yaml_options);

% profile clear
% profile -memory on -historysize 25000000



%% ---------------------------------------------------------------------- %
% Run sensitivity                                                         %
% ---------------------------------------------------------------------- %%

fn_sens(model_options);



%% ---------------------------------------------------------------------- %
% Run TFM                                                                 %
% ---------------------------------------------------------------------- %%

% scat_x = linspace(25.0e-3, 60.0e-3, 100);

% for scat_pos = 1:100
%     new_options.savename = sprintf('TFM SDH Scan %d', scat_pos);
%     new_options.scat_info = fn_scat_info( ...
%         "sdh", ...
%         1.0e-3, ...
%         6320/5e6, ...
%         3130/5e6, ...
%         deg2rad(0), ...
%         [[scat_x(scat_pos), 0.0, 27.5e-3]], ...
%         'ang_pts_over_2pi', 120 ...
%     );
% 
%     model_options = fn_default_model_options(new_options);
% 
%     fn_tfm(model_options);
% end
    

% new_options.savename = sprintf('TFM SDH Immersion Backwall Signal Testing');
% model_options = fn_default_model_options(new_options);
% 
% fn_tfm(model_options);



% if model_options.savepath ~= ""
%     cd(model_options.savepath)
%     if ~exist("profile", 'dir')
%         mkdir("profile")
%     end
%     profsave(profile('info'), "./profile")
% end
% profile report



% quit;
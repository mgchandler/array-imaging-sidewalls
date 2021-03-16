clear

new_options.probe_pitch =  1.00e-3;
new_options.pixel = 5.0e-3;
new_options.savepath = 'C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\array-imaging-sidewalls matlab\home-output\SDH Scan with weird geometry';
new_options.savename = 'SDH Refl - Valid Paths';
new_options.geom_shape.xmin = -25.0e-3;
new_options.geom_shape.xmax =  40.0e-3;
new_options.geom_shape.zmin =   0.0;
new_options.geom_shape.zmax =  45.0e-3;
% new_options.geometry = fn_make_geometry(0, 500, ...
% ...%     [new_options.geom_shape.xmax, 0.0, -1e-5], [new_options.geom_shape.xmin, 0.0, -1e-5], ...
%     [new_options.geom_shape.xmax, 0.0, new_options.geom_shape.zmax+1e-5], [new_options.geom_shape.xmin, 0.0, new_options.geom_shape.zmax+1e-5] ...
% );
new_options.geometry = fn_make_geometry(1, 500, ...
    [-25.00e-3, 0.0, 25.01e-3], [24.99e-3, 0.0, 25.01e-3], ...
    [ 24.99e-3, 0.0, 45.01e-3], [40.01e-3, 0.0, 45.01e-3], ...
    [ 40.01e-3, 0.0,  0.01e-3] ...
);
new_options.scat_info = fn_scat_info( ...
    "sdh", ...
    1.0e-3, ...
    6320/5e6, ...
    3130/5e6, ...
    deg2rad(0), ...
    [[32.5e-3, 0.0, 27.5e-3]], ...
    'ang_pts_over_2pi', 120 ...
);
% new_options.probe_standoff = 10e-3;
new_options.max_no_reflections = 1;
new_options.model_geometry = 0;
new_options.wall_for_imaging = "S2";
model_options = fn_default_model_options(new_options);

profile clear
profile -memory on



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
    

% new_options.savename = sprintf('TFM SDH Backwall Immersion');
% model_options = fn_default_model_options(new_options);

% fn_tfm(model_options);

profile report



% quit;
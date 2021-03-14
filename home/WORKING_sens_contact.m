clear

new_options.probe_pitch =  1.00e-3;
new_options.pixel = 5.0e-3;
new_options.savepath = 'C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\array-imaging-sidewalls matlab\home-output\Immersion Testing';
new_options.savename = sprintf('Sens SDH Backwall Immersion');
new_options.geom_shape.xmin = -25.0e-3;
new_options.geom_shape.xmax =  40.0e-3;
new_options.geom_shape.zmin =   0.0;
new_options.geom_shape.zmax =  45.0e-3;
% new_options.geometry = fn_make_geometry(0, 500, ...
% ...%     [new_options.geom_shape.xmax, 0.0, -1e-5], [new_options.geom_shape.xmin, 0.0, -1e-5], ...
%     [new_options.geom_shape.xmax, 0.0, new_options.geom_shape.zmax+1e-5], [new_options.geom_shape.xmin, 0.0, new_options.geom_shape.zmax+1e-5] ...
% );
new_options.geometry = fn_make_geometry(1, 500, ...
    [-25.00e-3, 0.0, 35.01e-3], [24.99e-3, 0.0, 35.01e-3], ...
    [ 24.99e-3, 0.0, 45.01e-3], [40.01e-3, 0.0, 45.01e-3], ...
    [ 40.01e-3, 0.0,  0.01e-3] ...
);
% new_options.probe_standoff = 10e-3;
new_options.max_no_reflections = 1;
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

new_options.savename = sprintf('TFM SDH Backwall Immersion');
model_options = fn_default_model_options(new_options);

% fn_tfm(model_options);

profile report



% quit;
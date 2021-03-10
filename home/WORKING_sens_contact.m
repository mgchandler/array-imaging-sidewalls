clear

model_config.PITCH =  1.00e-3;
model_config.PIXEL = 7.5e-3;
model_config.WALLS = 500;
model_config.VIEWS = 1;
model_config.GEOM = 0;
model_config.SETUP = 0;

new_options.savepath = 'C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\array-imaging-sidewalls matlab\home-output\Profiler Outputs\21.03.03 - Default setup';
new_options.savename = sprintf('Sens SDH Backwall for presentation');
new_options.geom_shape.xmin = -25.0e-3;
new_options.geom_shape.xmax =  25.0e-3;
new_options.geom_shape.zmin =   0.0;
new_options.geom_shape.zmax =  40.0e-3;
new_options.geometry = fn_make_geometry(0, 500, ...
    [new_options.geom_shape.xmax, 0.0, new_options.geom_shape.zmax+1e-5], [new_options.geom_shape.xmin, 0.0, new_options.geom_shape.zmax+1e-5] ...
);
new_options.max_no_reflections = 0;
model_options = fn_default_model_options(model_config, new_options);



%% ---------------------------------------------------------------------- %
% Run sensitivity                                                         %
% ---------------------------------------------------------------------- %%

profile clear
profile -memory on

fn_sens(model_config, model_options);

profile report



%% ---------------------------------------------------------------------- %
% Run TFM                                                                 %
% ---------------------------------------------------------------------- %%

% fn_tfm(model_config, model_options);



% quit;
model_config.PITCH =  1.00e-3;
model_config.PIXEL = 5.0e-3;
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
model_options = fn_default_model_options(new_options);

%% ---------------------------------------------------------------------- %
% Run sensitivity                                                         %
% ---------------------------------------------------------------------- %%

fn_sens(model_config, model_options);



%% ---------------------------------------------------------------------- %
% Run TFM                                                                 %
% ---------------------------------------------------------------------- %%

% fn_tfm(model_config, model_options);



% quit;
% Template script for running contact and immersion TFMs and sensitivities
% on a personal Matlab instance (i.e. not on BluePebble).

clear

%% ---------------------------------------------------------------------- %
% Setup input parameters                                                  %
% ---------------------------------------------------------------------- %%

model_config.PITCH =  1.00e-3;
model_config.PIXEL = 2.0e-3;
model_config.WALLS = 500;
model_config.VIEWS = 1;
model_config.GEOM = 0;
model_config.SETUP = 1;

new_options.scat_info = fn_scat_info("point", 6320.0, 3130.0, [[0.0, 0.0, 22.0e-3]]);
new_options.probe_standoff = 10e-3;
model_options = fn_default_model_options(new_options);



%% ---------------------------------------------------------------------- %
% Run sensitivity                                                         %
% ---------------------------------------------------------------------- %%

% fn_sens(model_config, model_options);



%% ---------------------------------------------------------------------- %
% Run TFM                                                                 %
% ---------------------------------------------------------------------- %%

fn_tfm(model_config, model_options);
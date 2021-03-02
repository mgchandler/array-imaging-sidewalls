function WORKING_high_res_sdh_back_side()
% Script which takes the array index from BluePebble as an input, works out
% exactly which geometry block this relates to, and performs the
% corresponding sidewall sensitivity calculation.
%
% INPUTS:
% - array_idx : integer
%       The array index PBS_ARRAY_INDEX in the shell script submitted to
%       BluePebble.
% - total_jobs : integer
%       The total number of array jobs submitted to BluePebble. Note that
%       there is no way of obtaining this number as a variable in the shell
%       script, so it must be manually coded into the shell.

% Move into engine directory.
% cd('/home/mc16535/Matlab_Sidewalls/array-imaging-sidewalls/engine')

% array_idx = 2;
% total_jobs = 100;
% VIEWS = 2;

% Assume that the submission directory is /array-imaging-sidewalls/bp
cd('../engine')

%% ---------------------------------------------------------------------- %
% Setup input parameters                                                  %
% ---------------------------------------------------------------------- %%

model_config.PITCH =  1.00e-3;
model_config.PIXEL = 1.0e-3;
model_config.WALLS = 1000;
model_config.VIEWS = 4;
model_config.GEOM = 0;
model_config.SETUP = 0;

% half_probe_width = 16*model_config.PITCH;
% full_xwidth = 100.0e-3;

% xmins = linspace(-half_probe_width, -full_xwidth, total_jobs);
% xmaxs = linspace(full_xwidth, half_probe_width, total_jobs);

new_options.scat_info = fn_scat_info( ...
    "sdh", ...
    1.0e-3, ...
    6320.0/5.0e+6, ...
    3130.0/5.0e+6, ...
    0, ...
    [[0, 0, 0]], ...
    'ang_pts_over_2pi', 120 ...
);
new_options.savepath = '/work/mc16535/Matlab_Sidewalls/array-sidewalls-imaging/';
% new_options.savepath = 'C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\array-imaging-sidewalls matlab\array-imaging-sidewalls\output';
new_options.savename = 'Back and Sidewall SDH';
new_options.geom_shape.xmin = -35.0e-3;
new_options.geom_shape.xmax =  35.0e-3;
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

quit;
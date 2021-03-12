function WORKING_valid_path_direct_sdh()
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
cd('../../engine')

%% ---------------------------------------------------------------------- %
% Setup input parameters                                                  %
% ---------------------------------------------------------------------- %%

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
new_options.savepath = '/work/mc16535/Matlab_Sidewalls/array-sidewalls-imaging/valid-paths';
% new_options.savepath = 'C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\array-imaging-sidewalls matlab\array-imaging-sidewalls\output';
new_options.savename = 'SDH Direct - Valid Paths';
new_options.geom_shape.xmin = -25.0e-3;
new_options.geom_shape.xmax =  40.0e-3;
new_options.geom_shape.zmin =   0.0;
new_options.geom_shape.zmax =  45.0e-3;
new_options.geometry = fn_make_geometry(1, 500, ...
    [-25.00e-3, 0.0, 25.01e-3], [24.99e-3, 0.0, 25.01e-3], ...
    [ 24.99e-3, 0.0, 45.01e-3], [40.00e-4, 0.0, 45.01e-3] ...
);
new_options.max_no_reflections = 0;
new_options.pixel = 1.0e-3;
model_options = fn_default_model_options(new_options);



%% ---------------------------------------------------------------------- %
% Run sensitivity                                                         %
% ---------------------------------------------------------------------- %%

fn_sens(model_options);



%% ---------------------------------------------------------------------- %
% Run TFM                                                                 %
% ---------------------------------------------------------------------- %%

% fn_tfm(model_config, model_options);

quit;
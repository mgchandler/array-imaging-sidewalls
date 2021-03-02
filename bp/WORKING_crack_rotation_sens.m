function WORKING_crack_rotation_sens(array_idx, total_jobs, VIEWS)
% Script which takes the array index from BluePebble as an input, works out
% exactly which rotation of crack this corresponds to, and runs the
% corresponding sensitivity code.
%
% INPUTS:
% - array_idx : integer
%       The array index PBS_ARRAY_INDEX in the shell script submitted to
%       BluePebble.
% - total_jobs : integer
%       The total number of array jobs submitted to BluePebble. Note that
%       there is no way of obtaining this number as a variable in the shell
%       script, so it must be manually coded into the shell.
% - VIEWS : integer
%       The value passed into the fn_sens function to determine which views
%       to image.

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

crack_angle = linspace(-pi, pi, total_jobs+1);
crack_angle = crack_angle(2:end);

model_config.PITCH =  1.00e-3;
model_config.PIXEL = 1.0e-3;
model_config.WALLS = 750;
model_config.VIEWS = VIEWS;
model_config.GEOM = 0;
model_config.SETUP = 0;

new_options.scat_info = fn_scat_info( ...
    "crack", ...
    6320.0, ...
    3130.0, ...
    2700.0, ...
    5.0e+6, ...
    2.0e-3, ...
    crack_angle(array_idx), ...
    [[0, 0, 0]], ...
    'nodes_per_wavelength', 20, ...
    'ang_pts_over_2pi', 120 ...
);
new_options.savepath = '/work/mc16535/Matlab_Sidewalls/array-sidewalls-imaging/crack-rotation';
% new_options.savepath = 'C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\array-imaging-sidewalls matlab\array-imaging-sidewalls\output';
new_options.savename = sprintf('Sens Contact %d Crack Angle %.2f', VIEWS, crack_angle(array_idx));
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
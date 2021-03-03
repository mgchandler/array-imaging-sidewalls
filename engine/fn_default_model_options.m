function model_options = fn_default_model_options(varargin)
% Returns the default model options if no inputs are given, and overwrites
% the defaults when these are provided.
%
% INPUTS:
% - model_options : OPTIONAL struct (1, 1)
%       Structure containing options which will overwrite the defaults.
%       Possible options to include as fields are:
%       - probe_els : integer : DEFAULT = 32
%       - probe_angle : double : DEFAULT = 0
%       - probe_standoff : double : DEFAULT = 0
%       - probe_frequency : double : DEFAULT = 5.0e+6
%       - el_length : double : DEFAULT = 0.0
%       - geom_shape : struct (1, 1) : DEFAULT = (-25e-3, 25e-3, 0e-3, 40e-3)
%       - multi_freq : logical : DEFAULT = 0
%       - material_params : struct : DEFAULT air and alum properties
%       - scat_info : struct : DEFAULT sdh located at (0, 22e-3)
%       - boxsize : integer : DEFAULT = 0
%       - savename : string : DEFAULT = "sens MODE GEOM VIEWS PITCH PIXEL WALLS"
%
% OUTPUTS:
% - model_options : struct (1, 1)
%       Structure containing the specified options.

model_options.probe_els = 32;
model_options.probe_angle = 0;
model_options.probe_standoff = 0;
model_options.probe_frequency = 5.0e+6;
model_options.el_length = 0;
model_options.geom_shape.xmin = -25.0e-3;
model_options.geom_shape.xmax =  25.0e-3;
model_options.geom_shape.zmin =   0.0e-3;
model_options.geom_shape.zmax =  40.0e-3;
model_options.multi_freq = 0;
model_options.material_params.couplant_speed = 340.0;
model_options.material_params.couplant_density = 1.2;
model_options.material_params.solid_long_speed = 6320.0;
model_options.material_params.solid_shear_speed = 3130.0;
model_options.material_params.solid_density = 2700;
model_options.scat_info = fn_scat_info( ...
    "sdh", ...
    1.0e-3, ...
    model_options.material_params.solid_long_speed/model_options.probe_frequency, ...
    model_options.material_params.solid_shear_speed/model_options.probe_frequency, ...
    deg2rad(0), ...
    [[0.0, 0.0, 22.0]], ...
    'ang_pts_over_2pi', 120 ...
);
model_options.boxsize = 0;
model_options.savepath = 'C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\array-imaging-sidewalls matlab\array-imaging-sidewalls\output';
model_options.savename = 'TFM-Sens Image Plot';

if nargin == 1
    new_options = varargin{1};
    option_fieldnames = fieldnames(new_options);
    for arg = 1:size(option_fieldnames, 1)
        field_val = getfield(new_options, option_fieldnames{arg});
        model_options = setfield(model_options, option_fieldnames{arg}, field_val);
    end
elseif nargin ~= 0
    error('fn_default_model_options: wrong number of inputs.')
end

end
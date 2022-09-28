function [Ims, Views_im, Views] = fn_tfm(model_options, varargin)
% Computes the sensitivity maps for an array in contact with the solid
% block being inspected. Currently works for a rectangular block with some
% depth defined as the zmax location, and sidewall location defined as the
% xmax location.
%
% INPUTS:
% - model_options : struct 
%       Structure containing options which will overwrite the defaults.
%       Possible options to include as fields are:
%       - data : struct : DEFAULT none
%           If FMC data generated externally, it can be passed in using
%           this option. The simulation will be skipped and the program
%           will proceed straight to the TFM.
%       - material : struct
%           - couplant_density : double : DEFAULT 1.2
%           - couplant_v : double : DEFAULT 340.0
%           - density : double : DEFAULT 2700.0
%           - modulus : double : DEFAULT 70e9
%           - poisson : double : DEFAULT 0.34
%       - mesh : struct
%           - geom : struct
%               - x : double array
%                   List of x coordinates
%               - y : double array
%                   List of z coordinates
%               - geometry : struct
%                   Output from fn_make_geometry()
%           - scat : struct
%               Output from fn_scat_info()
%       - model : struct
%           - boxsize : double : DEFAULT 0.0
%           - db_range : double : DEFAULT 40.0
%           - max_no_reflections : integer : DEFAULT 1
%           - model_geom : logical : DEFAULT 1
%           - multi_freq : logical : DEFAULT 0
%           - norm_to : double : DEFAULT 0
%           - npw : double : DEFAULT 45
%           - pixel : double : DEFAULT 0.5e-3
%           - savename : string : DEFAULT "TFM-Sens Image Plot"
%           - savepath : string : DEFAULT ""
%           - wall_for_imaging : string : DEFAULT "B1"
%       - probe : struct
%           - angle : double : DEFAULT 0
%           - freq : double : DEFAULT 5.0e+6
%           - num_els : integer : DEFAULT 32
%           - standoff : double : DEFAULT 0
%           - separation : double : DEFAULT 0.05e-3
%           - width : double : DEFAULT 0.45e-3

%% ---------------------------------------------------------------------- %
% Unpack model_config and model_options                                   %
% ---------------------------------------------------------------------- %%

savename = model_options.model.savename;
geometry = model_options.mesh.geom.geometry;

probe_angle = model_options.probe.angle;
probe_standoff = model_options.probe.standoff;
no_cycles = model_options.probe.cycles;
frequency = model_options.probe.freq;

no_walls = size(geometry, 1);

% Check whether we are in contact or immersion. If we are in contact, there
% will be no frontwall, and probe_standoff and probe_angle must equal zero.
% If we are in immersion, there must be a frontwall and probe_standoff must
% be non-zero.
is_frontwall = 0;
for wall = 1:no_walls
    if geometry(wall).name == "F"
        is_frontwall = 1;
        break
    end
end
if and(is_frontwall, probe_standoff ~= 0)
    is_contact = false;
elseif and(~is_frontwall, and(probe_standoff == 0, probe_angle == 0))
    is_contact = true;
else
    error('fn_sens: Invalid setup.')
end



%% ---------------------------------------------------------------------- %
% Simulate FMC data                                                       %
% ---------------------------------------------------------------------- %%
% If no FMC data supplied, then skip simulation.
if ~isstruct(model_options.data)

    [FMC_time, FMC_time_data] = fn_simulate_fmc(model_options);


%% ---------------------------------------------------------------------- %
% Load FMC data                                                           %
% ---------------------------------------------------------------------- %%
else
    
    FMC_time = model_options.data.time;
    FMC_time_data = model_options.data.data;
    
end



%% ---------------------------------------------------------------------- %
% Imaging                                                                 %
% ---------------------------------------------------------------------- %%

if nargin == 1
    [Ims, Views_im, Views] = fn_image_tfm(FMC_time, FMC_time_data, model_options);
else
    Views_im = varargin{1};
    [Ims, Views_im, Views] = fn_image_tfm(FMC_time, FMC_time_data, model_options, Views_im);
end



end
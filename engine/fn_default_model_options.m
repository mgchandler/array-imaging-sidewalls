function model_options = fn_default_model_options(varargin)
% Returns the default model options if no inputs are given, and overwrites
% the defaults when these are provided. Note that the architecture is laid
% out such that it is compatible with the AbaqusInputFileGeneration .yaml
% files.
%
% INPUTS:
% - new_options : OPTIONAL struct (1, 1)
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
%           - sdh : struct
%               - info : struct
%                   Output from fn_scat_info()
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
%
% OUTPUTS:
% - model_options : struct (1, 1)
%       Structure containing the specified options.

model_options.data = 0;

model_options.material.couplant_density = 1.2;
model_options.material.couplant_v = 340.0;
model_options.material.density = 2700.0;
model_options.material.modulus = 70.0e+9;
model_options.material.poisson = 0.34;

model_options.mesh.geom.shape.xmax =  25.0e-3;
model_options.mesh.geom.shape.xmin = -25.0e-3;
model_options.mesh.geom.shape.zmax =  40.0e-3;
model_options.mesh.geom.shape.zmin =   0.0e-3;
model_options.mesh.geom.geometry = fn_make_geometry(1, 1000, ...
        [model_options.mesh.geom.shape.xmin, 0.0, model_options.mesh.geom.shape.zmax], ...
        [model_options.mesh.geom.shape.xmax, 0.0, model_options.mesh.geom.shape.zmax], ...
        [model_options.mesh.geom.shape.xmax, 0.0, model_options.mesh.geom.shape.zmin] ...
);

model_options.model.boxsize = 0.0;
model_options.model.db_range = 40.0;
model_options.model.max_no_reflections = 1;
model_options.model.model_geom = 1;
model_options.model.multi_freq = 0;
model_options.model.norm_to = 0;
model_options.model.npw = 0;
model_options.model.pixel = 0.5e-3;
model_options.model.savename = "TFM-Sens Image Plot";
model_options.model.savepath = "";
model_options.model.wall_for_imaging = "B1";

model_options.probe.angle = 0;
model_options.probe.freq = 5.0e+6;
model_options.probe.num_els = 32;
model_options.probe.standoff = 0;
model_options.probe.separation = 0.05e-3;
model_options.probe.width = 0.45e-3;

model_options.mesh.sdh.info = fn_scat_info( ...
        "sdh", ...
        1.0e-3, ...
        6317.01/model_options.probe.freq, ...
        3110.28/model_options.probe.freq, ...
        deg2rad(0), ...
        [[10.0e-3, 0.0, 22.0e-3]], ...
        'ang_pts_over_2pi', 120 ...
);



if nargin == 1
    new_options = varargin{1};
    option_fieldnames = fieldnames(new_options);
    for arg = 1:size(option_fieldnames, 1)
        if ~isfield(model_options, option_fieldnames{arg})
            model_options.(option_fieldnames{arg}) = new_options.(option_fieldnames{arg});
        elseif ~isstruct(model_options.(option_fieldnames{arg}))
            model_options.(option_fieldnames{arg}) = new_options.(option_fieldnames{arg});
        else
            model_fieldnames = fieldnames(new_options.(option_fieldnames{arg}));
            for arg1 = 1:size(model_fieldnames, 1)
                model_options.(option_fieldnames{arg}).(model_fieldnames{arg1}) = new_options.(option_fieldnames{arg}).(model_fieldnames{arg1});
            end
        end
    end
elseif nargin ~= 0
    error('fn_default_model_options: wrong number of inputs.')
end

if isfield(model_options.mesh, 'geom')
    geom_corners = [cell2mat(model_options.mesh.geom.x); zeros(size(model_options.mesh.geom.x)); cell2mat(model_options.mesh.geom.y)].';
    model_options.mesh.geom.geometry = fn_make_geometry(1, 1000, ...
        geom_corners ...
    );
end

if model_options.model.max_no_reflections > size(model_options.mesh.geom.geometry, 1)
    error('fn_default_model_options: too many reflections for walls provided')
end

if isstruct(model_options.data)
    if ~and(isfield(model_options.data, 'time'), isfield(model_options.data, 'data'))
        error('fn_default_model_options: FMC_data field requires both time and data fields.')
    elseif size(model_options.data.time, 1) ~= size(model_options.data.data, 1)
        error('fn_default_model_options: FMC_data field requires both time and data fields be same size.')
    end
elseif model_options.data ~= 0
    error('fn_default_model_options: FMC_data must be zero or struct.')
end

end
function model_options = fn_default_model_options(varargin)
% Returns the default model options if no inputs are given, and overwrites
% the defaults when these are provided.
%
% INPUTS:
% - new_options : OPTIONAL struct (1, 1)
%       Structure containing options which will overwrite the defaults.
%       Possible options to include as fields are:
%       - pixel : double : DEFAULT = 5.0e-3
%           The pixel resolution of the resulting image.
%       - probe_els : integer : DEFAULT = 32
%           The number of elements in the probe. Note that the probe is
%           centered at the coordinates (0, 0, 0) if there is no standoff.
%       - probe_angle : double : DEFAULT = 0
%           The angle the probe makes with the front wall in radians. Note 
%           that in the contact case, this must be set to zero. In the 
%           immersion case, the standoff must be at least large enough to
%           account for the swing angle made by the probe when angle is
%           non-zero.
%       - probe_standoff : double : DEFAULT = 0
%           The standoff of the probe with respect to the front wall. Note
%           that in the contact case, this must be set to zero. In the 
%           mmersion case, the standoff must be at least large enough to
%           account foe the swing angle made by the probe when angle is
%           non-zero.
%       - probe_frequency : double : DEFAULT = 5.0e+6
%           The centre frequency of the probe.
%       - probe_pitch : double : DEFAULT = 1.0e-3
%           The pitch of the probe array. Must be positive real.
%       - el_length : double : DEFAULT = 0.0
%           The additional length each probe element on top of the probe
%           pitch.
%       - geom_shape : struct (1, 1)
%           The rectangular area which will be imaged at the end. This must
%           fully contain the geometry of the inspection block. Fields are
%           - xmin : double : DEFAULT = -25.0e-3
%           - xmax : double : DEFAULT =  25.0e-3
%           - zmin : double : DEFAULT =   0.0e-3
%           - zmax : double : DEFAULT =  40.0e-3
%       - multi_freq : logical : DEFAULT = 0
%           Logical switch for whether to use the multi-frequency model or
%           not. Default behaviour is to use the single-frequency model.
%           Contact requires the single-frequency model.
%       - material_params : struct (1, 1)
%           A structure detailing the various properties of the couplant
%           and the inspection block. Fields are
%           - couplant_speed : double : DEFAULT = 340.0
%           - couplant_density : double : DEFAULT = 1.2
%           - solid_long_speed : double : DEFAULT = 6320.0
%           - solid_shear_speed : double : DEFAULT = 3130.0
%           - solid_density : double : DEFAULT = 2700.0
%       - scat_info : struct : DEFAULT sdh located at (0, 22e-3)
%           Details on the scatterer which will be modelled. Note that
%           depending on the type of scatterer, different fields are
%           required. Should be an output from the fn_scat_info function.
%       - boxsize : double : DEFAULT = 0
%           The size of the box in metres which is drawn around the defect
%           in the TFM images. Has no effect in sensitivity.
%       - savepath : string : DEFAULT = ""
%           The path where the .fig and .mat files will be saved to. Note
%           that if equal to "", then no files will be saved. The
%           `savename` field will be unused.
%       - savename : string : DEFAULT = "sens MODE GEOM VIEWS PITCH PIXEL WALLS"
%           The name which will files will be saved as. The file types
%           (.fig and .mat) will be appended to this string.
%       - max_no_reflections : integer : DEFAULT = 1
%           The maximum number of reflections which will be made from the
%           walls of the geometry when tracing rays from probe to
%           scatterer.
%       - model_geometry : logical : DEFAULT = 0
%           Logical switch for whether to include the signal reflected from
%           the front and backwall. Note that this is not checked for in
%           fn_sens, as this signal is never modelled, and is only used in
%           fn_tfm.
%       - geometry : struct (no_walls, 1) : DEFAULT backwall
%           A structure containing all of the walls in the geometry, which
%           will be passed into ray tracing steps. Note that the presence
%           of a frontwall in this object defines whether we are in the
%           contact or immersion case.
%       - wall_for_imaging : string : DEFAULT 'B1'
%           When max_no_reflections > 1 and there is more than one wall,
%           imaged views will include reflections from the wall with name
%           equal to this string value.
%       - norm_to : double : DEFAULT 0
%           The reference value which will be used when normalising the
%           images. A value of zero is interpretted to mean that the
%           absolute maximum intensity across all views will be used;
%           otherwise the dB scale will be with reference to this value
%           (i.e. 0dB => intensity = norm_to).
%       - db_range_for_output : DEFAULT 40
%           The range over which the TFMs will be plotted (i.e. [-db_range,
%           0].
%       - FMC_data : struct : DEFAULT none
%           If FMC data generated externally, it can be passed in using
%           this option. The simulation will be skipped and the program
%           will proceed straight to the TFM.
%
% OUTPUTS:
% - model_options : struct (1, 1)
%       Structure containing the specified options.

model_options.pixel = 5.0e-3;
model_options.probe_els = 32;
model_options.probe_angle = 0;
model_options.probe_standoff = 0;
model_options.probe_frequency = 5.0e+6;
model_options.probe_pitch = 1.0e-3;
model_options.el_length = 0;

xmin = -25.0e-3;
xmax =  25.0e-3;
zmin =   0.0e-3;
zmax =  40.0e-3;

model_options.geom_shape.xmin = xmin;
model_options.geom_shape.xmax = xmax;
model_options.geom_shape.zmin = zmin;
model_options.geom_shape.zmax = zmax;
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
    [[10.0e-3, 0.0, 22.0e-3]], ...
    'ang_pts_over_2pi', 120 ...
);
model_options.boxsize = 0;
model_options.savepath = "";
model_options.savename = 'TFM-Sens Image Plot';
model_options.max_no_reflections = 0;
model_options.model_geometry = 0;
model_options.wall_for_imaging = 'B1';
model_options.norm_to = 0;
model_options.db_range_for_output = 40;

model_options.FMC_data = 0;

if nargin == 1
    new_options = varargin{1};
    option_fieldnames = fieldnames(new_options);
    for arg = 1:size(option_fieldnames, 1)
        model_options.(option_fieldnames{arg}) = new_options.(option_fieldnames{arg});
    end
elseif nargin ~= 0
    error('fn_default_model_options: wrong number of inputs.')
end

if ~isfield(model_options, 'geometry')
    model_options.geometry = fn_make_geometry(1, 500, ...
        [xmin, 0.0, zmax], [xmax, 0.0, zmax], [xmax, 0.0, zmin] ...
    );
end

if model_options.max_no_reflections > size(model_options.geometry, 1)
    error('fn_default_model_options: too many reflections for walls provided')
end

if isstruct(model_options.FMC_data)
    if ~and(isfield(model_options.FMC_data, 'time'), isfield(model_options.FMC_data, 'data'))
        error('fn_default_model_options: FMC_data field requires both time and data fields.')
    elseif size(model_options.FMC_data.time, 1) ~= size(model_options.FMC_data.data, 1)
        error('fn_default_model_options: FMC_data field requires both time and data fields be same size.')
    end
elseif model_options.FMC_data ~= 0
    error('fn_default_model_options: FMC_data must be zero or struct.')
end

if model_options.el_length > model_options.probe_pitch
    error('fn_default_model_options: Element length must not exceed element pitch.')
end

end
function path_info = fn_path_info(name, modes, geometries, speeds, mat_speeds, walls, medium_ids, densities, probe_frequency, pitch, probe_coords)
% Collects the information passed in into a single object which can then be
% passed into the ray computation function.
%
% INPUTS:
% - name : string
%       The name of the path.
% - modes : array (no_legs, 1)
%       Array of logicals. Each logical is a test of whether the
%       corresponding leg is shear or not (1 => shear, 0 => long).
% - geometries : array (no_walls, wall_pts, 3)
%       Contains the discretised 3D coordinates of each wall, provided in
%       the order the ray passes through them. Note that this variable can
%       alternatively be set to 0 if this is a direct contact path.
% - speeds : array (no_legs, 1)
%       The speed of the ray for each leg.
% - mat_speeds : array (3, 1)
%       Array of three speeds. In order, must be [couplant_speed,
%       solid_long_speed, solid_shear_speed]. Note that couplant speed must
%       be supplied even in the contact case.
% - walls : array (no_walls, 1)
%       Contains integers which identify the corresponding wall through
%       which the ray passes. 0 means no wall (direct contact), 1 means
%       frontwall, 2 means backwall, 3 means sidewall. Note that while they
%       should not be provided here, -1 means probe, -2 means scatterer.
% - medium_ids : array (no_legs, 1)
%       Contains logicals for each leg of the ray, treated as a logical
%       test of whether the leg is passing through a solid or not
% - densities : array (2, 1)
%       Densities of the couplant and the solid. Note that the couplant
%       must be supplied even in the contact case.
% - probe_frequency : double
%       Frequency of the probe inpute signal. Must be a single frequency:
%       frequencies for the multi-frequency model are supplied to the
%       fn_compute_ray function.
% - pitch : double
%       Pitch of the array.
% - probe_coords : array (probe_els, 3)
%       3D cooridnates of the probe elements.
%
% OUTPUTS:
% path_info : struct (1, 1)
%       Simple structure which just contains all of the supplied
%       information. This reduces the number of inputs to all subsequent
%       functions which require all of this information!

path_info.name = name;
path_info.modes = modes;
path_info.geometries = geometries;
path_info.speeds = speeds;
path_info.material_speeds = mat_speeds;
path_info.walls = walls;
path_info.medium_ids = medium_ids;
path_info.densities = densities;
path_info.probe_frequency = probe_frequency;
path_info.probe_pitch = pitch;
path_info.probe_coords = probe_coords;

end
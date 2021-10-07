function TR = fn_TR_coeff(medium_ids, mode_ids, inc_angles, material_speeds, densities)
% Computes the overall transmission/reflection coefficient for a ray as it
% travels from probe to scatterer. Only computes for a single ray, so this
% function assumes that the loop over probe and scatterer is in the parent
% function.
%
% INPUTS:
% - wall_ids : array (no_walls, 1)
%       Index conveying information on whether each wall the ray interfaces
%       with is a frontwall, backwall or sidewall. This is used to
%       determine how to compute the incident and outgoing angles.
% - medium_ids : array (no_legs, 1)
%       Identity of the media of each leg of the ray. Each element can be
%       treated as a logical test of whether the current medium is solid.
% - mode_ids : array (no_legs, 1)
%       Identity of the mode of each leg of the ray. Each element can be
%       treated as a logical test of whether the current mode is
%       transverse.
% - inc_angles : array (no_legs, 1)
%       Incident angles on each wall for each leg of the ray (note that the
%       last leg will not have an incident angle, so this one is omitted).
% - material_speeds : array (3, 1)
%       Array of speeds at which the ray travels in each medium. In order,
%       the array must contain the speed in the liquid; the longitudinal
%       speed in the solid; and the transverse speed in the solid.
% - densities : array (2, 1)
%       Array of densities of the media through which the ray travels. In
%       order, array must contain the density of the liquid; and the
%       density of the solid.
%
% OUTPUTS:
% - TR : complex
%       The total calculated transmission/reflection coefficient.

no_walls = size(inc_angles, 1);

TR = 1.0;

for wall = 1:no_walls
    if medium_ids(wall) == medium_ids(wall+1) % If the current medium is the
        % same as the next medium, then we must have a reflection.
        is_transmission = 0;
    else % The current medium is not the same, so we have a transmission.
        is_transmission = 1;
    end
        
    TR = TR * fn_transrefl( ...
        is_transmission, medium_ids(wall), mode_ids(wall), mode_ids(wall+1), ...
        inc_angles(wall), material_speeds, densities ...
    );
end

end
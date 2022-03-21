function directivity = fn_line_directivity(theta, long_lambda, shear_lambda, mode, c44)
% Computes the inherent directivity of the probe in contact with a solid,
% assuming an infinitely long strip source with finite width. Functions
% taken from Miller & Pursey (1954): Eqs 74, 93 and 94.
%
% INPUTS:
% - theta : double
%       angle of the ray with respect to the normal from the probe
%       axis.
% - long_lambda : double
%       wavelength of the longitudinal wave in the solid.
% - shear_lambda : double
%       wavelength of the shear wave in the solid.
% - mode : logical
%       mode of the ray for which directivity will be computed. Treated
%       as a logical test of whether ray is shear.
% - c44 : double
%       shear modulus of the material.
%
% OUTPUTS:
% - directivity : complex
%       The line directivity calculated for the provided mode.
%
% REFERENCES:
% [1] : G. Miller, H. Pursey, et al., The field and radiation impedance of
%       mechanical radiators on the free surface of a semi-infinite
%       isotropic solid, in Proceedings of the Royal Society of London,
%       Series A, Mathematical and Physical Sciences, vol. 223, 521-541,
%       (1954), doi:10.1098/rspa.1954.0134.

mu = long_lambda / shear_lambda;
sin_theta = sin(theta);
cos_theta = cos(theta);

if ~mode % If we want the longitudinal directivity, Eq 93.
    directivity = ( ...exp(3i/4 * pi) * sqrt(2/pi) / c44 * ( ...
        cos_theta * (mu^2 - 2 * sin_theta^2) / F_0(sin_theta, mu) ...
    );
else % Then we want the shear directivity, Eq 94.
    directivity = ( ...exp(5i/4 * pi) * sqrt(2/pi) / c44 * ( ...
        mu^2.5 * ...
        sin(2*theta) * sqrt(mu^2 * sin_theta^2 - 1) / F_0(sin_theta*mu, mu) ...
    );
end

end

function q = F_0(z, mu) % Helper function, Eq 74.
q = (2 * z^2 - mu^2)^2 - 4 * z^2 * sqrt(z^2 - 1) * sqrt(z^2 - mu^2);
end
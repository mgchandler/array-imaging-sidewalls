function beam_spread = fn_beamspread_2d(dists, alphas, speeds, varargin)
% Computes the beam spreading in 2D by finding the virtual distance of the
% ray. The function assumes that the ray has already been found. Note that
% this function was designed to only find one beam spread value: loops
% around probe element, scatterer and frequency should occur in the parent
% funciton if required.
%
% INPUTS:
% - dists : array (no_legs, 1)
%       The Euclidean distance of each leg in the array. When calculating
%       the forward beamspread (i.e. probe -> scatterer), the first element
%       should be the distance of the first leg (starts at the probe), and
%       the last element should be the last leg (ends at the scatterer).
% - alphas : array (no_walls, 1)
%       Incident angle that each leg of the ray makes with the following
%       wall. Order must correspond to the order of `dists` array. Note
%       that this array is one element shorter than `dists`: the incident
%       angle on the scatterer is not required.
% - speeds : array (no_legs, 1)
%       Speed of the ray for each leg. Order must be the same as `dists`.
%
% OUTPUTS:
% - beam_spread : double
%       bs = 1 / sqrt(virtual_distance), method taken from the Appendix of
%       Budyn [1]. Virtual distance is the distance which the ray appears
%       to have travelled in the last leg of the ray, based on the virtual
%       projection of the beam (distance M2' -> M3 in Fig 7b in source).
%
% REFERENCES:
% [1] : N. Budyn, R. Bevan, et al., "A Model for Multiview Ultrasonic Array
%       Inspection of Small Two-DImensional Defects," in IEEE Transactions
%       on Ultrasonics, Ferroelectrics, and Frequency Control, vol. 66, no.
%       6, 1129-1139, Jun 2019 (doi:10.1109/TUFFC.2019.2909988)

if nargin > 0
    freq = varargin{1};
    last_k = 2 * pi * freq / speeds(end);
else
    last_k = 1;
end

[no_legs, ~] = size(dists);
gamma_list = zeros(no_legs-1, 1);

virtual_distance = dists(1);

% If there is more than one leg in the virtual distance.
if no_legs ~= 1
    for k = 1:no_legs-1
        nu = speeds(k) / speeds(k+1);
        sin_alpha = sin(alphas(k));
        cos_alpha = cos(alphas(k));

        gamma_list(k) = (nu*nu - sin_alpha*sin_alpha) / (nu*cos_alpha*cos_alpha);
    end

    for k = 1:no_legs-1
        r = dists(k+1);
        gamma = 1.0;
        for ii = 1:k
            gamma = gamma*gamma_list(ii);
        end
        virtual_distance = virtual_distance + r/gamma;
    end
end

% beam_spread = 1/sqrt(last_k * virtual_distance);
beam_spread = 1/sqrt(virtual_distance);

end
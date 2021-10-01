function scat_amp = fn_sdh_scatterer_amp(inc_angle, out_angle, sdh_radius, lambda_L, lambda_T, inc_mode, out_mode)
% Computes the scattering amplitude for a side drilled hole (SDH) for some
% incident and outgoing angle.
%
% INPUTS:
% - inc_angle : double OR array (1, ang_pts_over_2pi)
%       Incident angle on the scatterer.
% - out_angle : double OR array (1, ang_pts_over_2pi)
%       Outgoing angle from the scatterer.
% - sdh_radius : double
%       Radius of the SDH.
% - lambda_L : double OR array (1, num_freqs)
%       Wavelength of the longitudinal mode.
% - lambda_T : double OR array (1, num_freqs)
%       Wavelength of the transverse mode.
% - inc_mode : logical
%       Logical test of whether the incident mode is shear.
% - out_mode : logical
%       Logical test of whether the outgoing mode is shear.
%
% OUTPUTS:
% - scat_amp : complex
%       Scattering amplitude from a side-drilled hole for the provided
%       modes and angles. Note: need to update to allow matrix
%       calcualtions.

%% Compute terms used in all scattering amplitudes.
if length(inc_angle) == 1
    theta = out_angle - inc_angle + pi; % Add pi as NDE convention is incident 
%       direction is angle between -ve incident wavevector and +ve x-axis.
else
    theta = (inc_angle + 2*pi)';
end

alpha = 2 * pi * sdh_radius ./ lambda_L;
beta = 2 * pi * sdh_radius ./ lambda_T;

N = max([10, ceil(4*alpha), ceil(4*beta)]); % Brind and Achenbach recommended truncation for 10E-6 relative accuracy.
n = [0:N];

epsilon_n = [1, 2*ones(1, N)]; % (2 - δ_0n)

C_n_1_alpha = C_n_i(n, 1, alpha, beta);
C_n_1_beta = C_n_i(n, 1, beta, beta);
D_n_1_alpha = D_n_i(n, 1, alpha);
D_n_1_beta = D_n_i(n, 1, beta);

%% Compute scattering amplitudes.
if ~inc_mode % If incident mode is longitudinal.
    if ~out_mode % If outgoing mode is longitudinal.
        C_n_2_alpha = C_n_i(n, 2, alpha, beta);
        D_n_2_alpha = D_n_i(n, 2, alpha);
        
        a = 1i / (2*alpha);
        numer = C_n_2_alpha .* C_n_1_beta - D_n_2_alpha .* D_n_1_beta;
        denom = C_n_1_alpha .* C_n_1_beta - D_n_1_alpha .* D_n_1_beta;

        A_n_L = 1i / (2 * alpha) * ( ...
            1 + (C_n_2_alpha .* C_n_1_beta - D_n_2_alpha .* D_n_1_beta) ./ ...
                (C_n_1_alpha .* C_n_1_beta - D_n_1_alpha .* D_n_1_beta)...
        );
    
%       Note that in Lopez equations, factor of lambda/sqrt(r) is
%       beamspreading, so we omit it here.
        scat_amp = sqrt(1i) / pi * sum( ...
            ones(size(theta)) * (epsilon_n .* alpha .* A_n_L) .* cos(theta * n), 2 ...
        );
        
    else % If incident mode is longitudinal, outgoing mode is shear.
        B_n_L = 2*n / (pi * alpha) .* ( ...
            (n.^2 - beta^2/2 - 1) ./ ...
            (C_n_1_alpha .* C_n_1_beta - D_n_1_alpha .* D_n_1_beta) ...
        );
    
        scat_amp = sqrt(1i) / pi* sum( ...
            ones(size(theta)) * (epsilon_n .* beta .* B_n_L) .* sin(theta * n), 2 ...
        );
        
    end
else % If incident mode is shear.
    if ~out_mode % If outgoing mode is longitudinal.
        A_n_T = 2*n / (pi * beta) .* ( ...
            (n.^2 - beta^2/2 - 1) ./ ...
            (C_n_1_alpha .* C_n_1_beta - D_n_1_alpha .* D_n_1_beta) ...
        );
    
        scat_amp = sqrt(1i) / pi * sum( ...
            ones(size(theta)) * (epsilon_n .* alpha .* A_n_T) .* sin(theta * n), 2 ...
        );
        
    else % If incident mode is shear, outgoing mode is shear.
        C_n_2_beta = C_n_i(n, 2, beta, beta);
        D_n_2_beta = D_n_i(n, 2, beta);

        B_n_T = 1i / (2 * beta) * ( ...
            1 + (C_n_2_beta .* C_n_1_alpha - D_n_2_beta .* D_n_1_alpha) ./ ...
                (C_n_1_alpha .* C_n_1_beta - D_n_1_alpha .* D_n_1_beta) ...
        );
    
        scat_amp = sqrt(1i) / pi * sum( ...
            ones(size(theta)) * (epsilon_n .* beta .* B_n_T) .* cos(theta * n), 2 ...
        );
        
    end
end

if length(inc_angle) ~= 1
    amp = zeros(size(scat_amp, 1));
    for jj = 1:size(scat_amp, 1)
        amp(:, jj) = scat_amp;
        scat_amp = circshift(scat_amp, 1);
    end
    scat_amp = amp;
end

end



function C = C_n_i(n, i, x, beta)
C = (n.^2 + n - beta^2/2) .* besselh(n, i, x) - x * besselh(n-1, i, x);
end

function D = D_n_i(n, i, x)
D = (n.^2 + n) .* besselh(n, i, x) - x * n .* besselh(n-1, i, x);
end
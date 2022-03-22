function scat_matrix = fn_scattering_matrix(scat_info, ang_pts_over_2pi)
% Computes scattering matrices for a single type of scatterer over a range
% of points.
%
% INPUTS:
% - scat_info : struct (1, 1)
%       Contains information on all scatterers which are being modelled.
%       Scatterer coordinates are in field `image_block`. Must be an output
%       from the fn_scat_info function.
% - ang_pts_over_2pi : integer
%       Number of points over 2 pi over which scattering amplitudes will be
%       precomputed.
%
% OUTPUTS:
% - scat_matrix : array (ang_pts_over_2pi, ang_pts_over_2pi)
%       Complex matrix. The (i, j)th position gives the scattering
%       amplitude with inc_angle(j) and out_angle(i), where inc and out
%       angles are linspaces from -pi to pi discretised into
%       ang_pts_over_2pi points.
    
% Create incident and outgoing angles.
inc_angle = linspace(-pi, pi, ang_pts_over_2pi+1);
inc_angle = inc_angle(1:end-1);
out_angle = linspace(-pi, pi, ang_pts_over_2pi+1);
out_angle = out_angle(1:end-1);

% Initialise matrices.
scat_matrix(1).name = "LL";
scat_matrix(1).matr = zeros(ang_pts_over_2pi);
scat_matrix(2).name = "LT";
scat_matrix(2).matr = zeros(ang_pts_over_2pi);
scat_matrix(3).name = "TL";
scat_matrix(3).matr = zeros(ang_pts_over_2pi);
scat_matrix(4).name = "TT";
scat_matrix(4).matr = zeros(ang_pts_over_2pi);
    
if scat_info.type == "sdh"
    % Unpack scatterer info.
    sdh_radius = scat_info.r(1);
    lambda_L = scat_info.lambdaL;
    lambda_T = scat_info.lambdaT;
    
    % Iterate. Need to change the sdh_scat_amp function to calculate by
    % matrix multiplication.
    scat_matrix(1).matr = fn_sdh_scatterer_amp( ...
        inc_angle, out_angle, sdh_radius, lambda_L, lambda_T, 0, 0 ...
    );
    scat_matrix(2).matr = fn_sdh_scatterer_amp( ...
        inc_angle, out_angle, sdh_radius, lambda_L, lambda_T, 0, 1 ...
    );
    scat_matrix(3).matr = fn_sdh_scatterer_amp( ...
        inc_angle, out_angle, sdh_radius, lambda_L, lambda_T, 1, 0 ...
    );
    scat_matrix(4).matr = fn_sdh_scatterer_amp( ...
        inc_angle, out_angle, sdh_radius, lambda_L, lambda_T, 1, 1 ...
    );
elseif scat_info.type == "crack"
    % Unpack scatterer info.
    density = scat_info.dens;
    vel_L = scat_info.vL;
    vel_T = scat_info.vT;
    frequency = scat_info.freq;
    crack_length = scat_info.crack_length;
    nodes_per_wavelength = scat_info.nodes_per_wavelength;
    
    % Iterate. Need to change the sdh_scat_amp function to calculate by
    % matrix multiplication.
    scat_matrix(1).matr = fn_crack_scatterer_amp( ...
        inc_angle, out_angle, density, vel_L, vel_T, ...
        frequency, 0, 0, crack_length, nodes_per_wavelength ...
    );
    scat_matrix(2).matr = fn_crack_scatterer_amp( ...
        inc_angle, out_angle, density, vel_L, vel_T, ...
        frequency, 0, 1, crack_length, nodes_per_wavelength ...
    );
    scat_matrix(3).matr = fn_crack_scatterer_amp( ...
        inc_angle, out_angle, density, vel_L, vel_T, ...
        frequency, 1, 0, crack_length, nodes_per_wavelength ...
    );
    scat_matrix(4).matr = fn_crack_scatterer_amp( ...
        inc_angle, out_angle, density, vel_L, vel_T, ...
        frequency, 1, 1, crack_length, nodes_per_wavelength ...
    );
else
    disp('Invalid scatterer type.')
end
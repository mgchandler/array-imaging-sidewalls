function scat_amps = fn_scattering_amps(view, freq_array)
% Computes amplitude of a scatterer for a specified view. Note that the
% current implementation assumes that all scatterers are the same type.
%
% INPUTS:
% - view : struct (1, 1)
%       Single view for which scattering amplitudes will be calculated.
%       Must be an output from the fn_create_view function.
% - freq_array : array (no_freqs, 1)
%       Non-zero frequencies for which scattering amplitudes will be
%       calculated. Note that if the multi-frequency model is being used,
%       then this should be the non-zero values output from the
%       fn_create_input_signal function.
%
% OUTPUTS:
% - scat_amps : array (probe_els^2, no_scats, no_freqs)
%       Complex scattering amplitudes calculated for an incoming and
%       outgoing angle, at each frequency for all scatterers.

% Unpack scatterer and path info.
path1 = view.path_1;
path2 = view.path_2;
path1_info = path1.path_info;
path2_info = path2.path_info;
scat_info = path1.scat_info;

[probe_els, num_scatterers] = size(path1.min_times);
[num_freqs, ~] = size(freq_array);

scat_amps = zeros(probe_els^2, num_scatterers, num_freqs);

if isfield(scat_info, 'matrix')
    inc_mode = path1_info.modes(end);
    out_mode = path2_info.modes(end);
    scat_matrix = scat_info.matrix(inc_mode*2 + out_mode + 1).matr;
    theta = linspace(-pi, pi, size(scat_matrix, 1));
    % Loop over element pairs. Could vectorise over element as well as
    % scatterers, but out-of-memory error arises for high element number
    % and small pixel size.
    for el = 1 : probe_els^2
        inc_angles = path1.weights.inc_theta(view.probe_txrx(el, 1), :, end, 2) - scat_info.angle;
        out_angles = path2.weights.inc_theta(view.probe_txrx(el, 2), :, 1, 2) - scat_info.angle;
        scat_amps(el, :, :) = interp2(theta, theta, scat_matrix, inc_angles, out_angles);
    end
else
    el = 1;
    for tx = 1 : probe_els
        for rx = 1 : probe_els
            for scat = 1 : num_scatterers
                for freq = 1 : num_freqs

                    if scat_info.type == "point"
                        inc_mode = path1_info.modes(end);
                        out_mode = path2_info.modes(end);
                        v_L = path1_info.material_speeds(2);
                        v_T = path1_info.material_speeds(3);

                        scat_amps(el, scat, freq) = fn_point_scatterer_amp( ...
                            inc_mode, out_mode, v_L, v_T ...
                        );



                    elseif scat_info.type == "sdh"
                        inc_angle_on_scat = path1.weights.inc_theta(tx, scat, end, 2) - scat_info.angle;
                        out_angle_on_scat = path2.weights.inv_out_theta(rx, scat, 1, 2) - scat_info.angle;

                        inc_mode = path1_info.modes(end);
                        out_mode = path2_info.modes(end);
                        lambda_L = path1_info.material_speeds(2) / freq_array(freq);
                        lambda_T = path1_info.material_speeds(3) / freq_array(freq);

                        scat_amps(el, scat, freq) = fn_sdh_scatterer_amp( ...
                            inc_angle_on_scat, out_angle_on_scat, scat_info.r, ...
                            lambda_L, lambda_T, inc_mode, out_mode ...
                        );



                    elseif scat_info.type == "crack"
                        inc_angle_on_scat = path1.weights.inc_theta(tx, scat, end, 2) - scat_info.angle;
                        out_angle_on_scat = path2.weights.inv_out_theta(rx, scat, 1, 2) - scat_info.angle;

                        density = scat_info.dens;
                        vel_L = scat_info.vL;
                        vel_T = scat_info.vT;
                        frequency = freq_array(freq);
                        inc_mode = path1_info.modes(end);
                        out_mode = path2_info.modes(end);
                        crack_length = scat_info.crack_length;
                        nodes_per_wavelength = scat_info.nodes_per_wavelength;

                        scat_amps(el, scat, freq) = fn_crack_scatterer_amp( ...
                            inc_angle_on_scat, out_angle_on_scat, ...
                            density, vel_L, vel_T, frequency, inc_mode, out_mode, ...
                            crack_length, nodes_per_wavelength ...
                        );



                    elseif scat_info.type == "image"
                        % We are imaging, and therefore do not require any
                        % scattering amplitude
                        continue



                    else
                        disp("Invalid scatterer type.")
                    end
                end
            end

            el = el + 1;
        end
    end
end

% el = 1;
% for tx = 1 : probe_els
%     for rx = 1 : probe_els
%         for scat = 1 : num_scatterers
%             for freq = 1 : num_freqs
%                 
%                 if scat_info.type == "point"
%                     inc_mode = path1_info.modes(end);
%                     out_mode = path2_info.modes(end);
%                     v_L = path1_info.material_speeds(2);
%                     v_T = path1_info.material_speeds(3);
%                     
%                     scat_amps(el, scat, freq) = fn_point_scatterer_amp( ...
%                         inc_mode, out_mode, v_L, v_T ...
%                     );
%                 
%                 
%                 
%                 elseif scat_info.type == "sdh"
%                     inc_angle_on_scat = path1.weights.inc_theta(tx, scat, end, 2) - scat_info.angle;
%                     out_angle_on_scat = path2.weights.inv_out_theta(rx, scat, 1, 2) - scat_info.angle;
% 
%                     if isfield(scat_info, 'matrix') % If matrices have been precomputed.
%                         % Get the relevant matrix for this view.
%                         inc_mode = path1_info.modes(end);
%                         out_mode = path2_info.modes(end);
%                         scat_matrix = scat_info.matrix(inc_mode*2 + out_mode + 1).matr;
%                         
%                         scat_amps(el, scat, freq) = fn_scattering_bilinear_interp(scat_matrix, inc_angle_on_scat, out_angle_on_scat);
% %                         scat_amps(el, scat, freq) = interp2(T1, T2, scat_matrix, inc_angle_on_scat, out_angle_on_scat, 'linear');
% 
%                     else % Matrices have not been found, instead compute using function.
%                         inc_mode = path1_info.modes(end);
%                         out_mode = path2_info.modes(end);
%                         lambda_L = path1_info.material_speeds(2) / freq_array(freq);
%                         lambda_T = path1_info.material_speeds(3) / freq_array(freq);
% 
%                         scat_amps(el, scat, freq) = fn_sdh_scatterer_amp( ...
%                             inc_angle_on_scat, out_angle_on_scat, scat_info.r, ...
%                             lambda_L, lambda_T, inc_mode, out_mode ...
%                         );
%                     end
%                     
%                     
%                     
%                 elseif scat_info.type == "crack"
%                     inc_angle_on_scat = path1.weights.inc_theta(tx, scat, end, 2) - scat_info.angle;
%                     out_angle_on_scat = path2.weights.inv_out_theta(rx, scat, 1, 2) - scat_info.angle;
% 
%                     if isfield(scat_info, 'matrix') % If matrices have been precomputed.
%                         % Get the relevant matrix for this view.
%                         inc_mode = path1_info.modes(end);
%                         out_mode = path2_info.modes(end);
%                         scat_matrix = scat_info.matrix(inc_mode*2 + out_mode + 1).matr;
%                         scat_amps(el, scat, freq) = fn_scattering_bilinear_interp(scat_matrix, inc_angle_on_scat, out_angle_on_scat);
% 
%                     else % Matrices have not been found, instead compute using function.
%                         density = scat_info.dens;
%                         vel_L = scat_info.vL;
%                         vel_T = scat_info.vT;
%                         frequency = freq_array(freq);
%                         inc_mode = path1_info.modes(end);
%                         out_mode = path2_info.modes(end);
%                         crack_length = scat_info.crack_length;
%                         nodes_per_wavelength = scat_info.nodes_per_wavelength;
%                         
%                         scat_amps(el, scat, freq) = fn_crack_scatterer_amp( ...
%                             inc_angle_on_scat, out_angle_on_scat, ...
%                             density, vel_L, vel_T, frequency, inc_mode, out_mode, ...
%                             crack_length, nodes_per_wavelength ...
%                         );
%                     end
%                     
%                     
%                     
%                 elseif scat_info.type == "image"
%                     % We are imaging, and therefore do not require any
%                     % scattering amplitude
%                     continue
%                     
%                     
%                     
%                 else
%                     disp("Invalid scatterer type.")
%                 end
%             end
%         end
%         
%         el = el + 1;
%     end
% end
end
function coeff = fn_transrefl(is_transmission, material_inc, mode_inc, mode_out, angle_inc, mat_speeds, density)
% Helper function for the fn_TR_coeff function which does all of the work.
%
% INPUTS:
% - is_transmission : logical
%       Logical test for transmission - set to 1 if the path is transmitted
%       through the current interface, and 0 if the path is reflected from.
% - material_inc : logical
%       Can be treated as a logical test for whether the incident material
%       is solid: set to 0 if incident medium is liquid, and 1 if solid.
% - mode_inc : logical
%       Can be treated as a logical test for whether the incident mode is
%       transverse: set to 0 if incident mode is longitudinal, and 1 if
%       transverse. Note: if material_inc is 0 (i.e. liquid), then mode_inc
%       must be 0 (longitudinal).
% - mode_out : logical
%       Can again be treated as a logical test for whether the output mode
%       is transverse.
% - angle_inc : double
%       Incident angle at the interface. Other angles can be determined
%       from the incident material, mode and whether it's a transmission/
%       reflection.
% - mat_speeds = [liquid_speed, solid_long_speed, solid_trans_speed] : double
%       The longitudinal speed in the liquid, the longitudinal speed in the
%       solid and the transverse speed in the solid respectively.
% - density = [liquid_dens, solid_dens] : double
%       The density of the liquid and the density of the solid
%       respectively.
%
% OUTPUTS:
% - coeff : complex
%       The calculated coefficient.

liquid_speed = mat_speeds(1);
solid_long_speed = mat_speeds(2);
solid_trans_speed = mat_speeds(3);
liquid_dens = density(1);
solid_dens = density(2);

% As equations are different depending on whether the ray is transmitted or
% reflected, work out which one this is.
if is_transmission % Do transmission. Output material is opposite to material_inc.
    material_out = ~material_inc;
    if ~material_inc % If incident material is liquid. Output material must be solid.
        alpha_fluid = angle_inc;
        alpha_l = conj(asin(sin(alpha_fluid) * solid_long_speed / liquid_speed));
        alpha_t = conj(asin(sin(alpha_fluid) * solid_trans_speed / liquid_speed));
        
        ct_cl2 = (solid_trans_speed * solid_trans_speed) / (solid_long_speed * solid_long_speed);
        cos_2_alpha_t = cos(2 * alpha_t);
        N = ( ...
            ct_cl2 * sin(2 * alpha_l) .* sin(2 * alpha_t) ...
            + cos_2_alpha_t .* cos_2_alpha_t ...
            + liquid_dens * liquid_speed / (solid_dens * solid_long_speed) * cos(alpha_l) ./ cos(alpha_fluid) ...
        );
    
        if ~mode_out % If output mode is longitudinal.
            % trans_l.
            coeff = 2.0 * cos_2_alpha_t ./ N;
            
            % For conversion to displacement: z = 
            % (rho_inc * vel_inc)/(rho_out * c_out)
            z = (liquid_dens * liquid_speed) / (solid_dens * solid_long_speed);
            coeff = z * coeff;
            
        else % Output mode must be transverse.
            % trans_t.
            coeff = -2.0 * ct_cl2 * sin(2 * alpha_l) ./ N;
            
            % For conversion to displacement: z = 
            % (rho_inc * vel_inc)/(rho_out * c_out)
            z = (liquid_dens * liquid_speed) / (solid_dens * solid_trans_speed);
            coeff = z * coeff;
            
        end
        
    else % Incident material must be solid. Incident mode can be longitudinal or transverse.
        if ~mode_inc % If incident mode is longitudinal
            alpha_l = angle_inc;
            alpha_fluid = conj(asin(sin(alpha_l) * liquid_speed / solid_long_speed));
            alpha_t = conj(asin(sin(alpha_l) * solid_trans_speed / solid_long_speed));
            
            ct_cl2 = (solid_trans_speed * solid_trans_speed) / (solid_long_speed * solid_long_speed);
            cos_2_alpha_t = cos(2 * alpha_t);
            N = ( ...
                ct_cl2 * sin(2 * alpha_l) .* sin(2 * alpha_t) ...
                + cos_2_alpha_t .* cos_2_alpha_t ...
                + liquid_dens * liquid_speed / (solid_dens * solid_long_speed) * cos(alpha_l) ./ cos(alpha_fluid) ...
            );
        
            % trans.
            coeff = ( ...
                2 ...
                * liquid_dens ...
                * liquid_speed ...
                * cos(alpha_l) ...
                .* cos(2 * alpha_t) ...
                ./ (N * solid_dens * solid_long_speed .* cos(alpha_fluid)) ...
            );
            
            % For conversion to displacement: z = 
            % (rho_inc * vel_inc)/(rho_out * c_out)
            z = (solid_dens * solid_long_speed) / (liquid_dens * liquid_speed);
            coeff = z * coeff;
            
        else % Incident mode must be transverse.
            alpha_t = angle_inc;
            alpha_fluid = conj(asin(sin(alpha_t) * liquid_speed / solid_trans_speed));
            alpha_l = conj(asin(sin(alpha_t) * solid_long_speed / solid_trans_speed));
            
            ct_cl2 = (solid_trans_speed * solid_trans_speed) / (solid_long_speed * solid_long_speed);
            cos_2_alpha_t = cos(2 * alpha_t);
            N = ( ...
                ct_cl2 * sin(2 * alpha_l) .* sin(2 * alpha_t) ...
                + cos_2_alpha_t .* cos_2_alpha_t ...
                + liquid_dens * liquid_speed / (solid_dens * solid_long_speed) * cos(alpha_l) ./ cos(alpha_fluid) ...
            );
        
            % trans.
            coeff = ( ...
                2 ...
                * liquid_dens ...
                * liquid_speed ...
                * cos(alpha_l) ...
                * sin(2 * alpha_t) ...
                ./ (N * solid_dens * solid_long_speed .* cos(alpha_fluid)) ...
            );
            
            % For conversion to displacement: z = 
            % (rho_inc * vel_inc)/(rho_out * c_out)
            z = (solid_dens * solid_trans_speed) / (liquid_dens * liquid_speed);
            coeff = z * coeff;
            
        end
    end
    
else % Do reflection. Output material is same as material_inc.
    material_out = material_inc;
    if material_inc % If incident material is solid. Output material must be solid.
        if ~mode_inc % If incident mode is longitudinal.
            alpha_l = angle_inc;
            alpha_fluid = conj(asin(sin(alpha_l) * liquid_speed / solid_long_speed));
            alpha_t = conj(asin(sin(alpha_l) * solid_trans_speed / solid_long_speed));

            ct_cl2 = (solid_trans_speed * solid_trans_speed) / (solid_long_speed * solid_long_speed);
            cos_2_alpha_t = cos(2 * alpha_t);
            N = ( ...
                ct_cl2 * sin(2 * alpha_l) .* sin(2 * alpha_t) ...
                + cos_2_alpha_t .* cos_2_alpha_t ...
                + liquid_dens * liquid_speed / (solid_dens * solid_long_speed) * cos(alpha_l) ./ cos(alpha_fluid) ...
            );
        
            if ~mode_out % If output mode is longitudinal.
                % refl_l.
                coeff = ( ...
                    ct_cl2 * sin(2 * alpha_l) .* sin(2 * alpha_t) ...
                    - cos_2_alpha_t .* cos_2_alpha_t ...
                    + liquid_dens * liquid_speed / (solid_dens * solid_long_speed) * cos(alpha_l) ./ cos(alpha_fluid) ...
                ) ./ N;
            
                % For conversion to displacement: z = 
                % (rho_inc * vel_inc)/(rho_out * c_out). Note that when reflecting,
                % if inc_mode == out_mode, then z = 1, so no conversion required.
                
            else % Output mode must be transverse.
                % refl_t.
                coeff = (2 * ct_cl2 * sin(2 * alpha_l) .* cos(2 * alpha_t)) ./ N;
            
                % For conversion to displacement: z = 
                % (rho_inc * vel_inc)/(rho_out * c_out)
                z = (solid_long_speed) / (solid_trans_speed);
                coeff = z * coeff;
                
            end
            
        else % Incident mode must be transverse.
            alpha_t = angle_inc;
            alpha_fluid = conj(asin(sin(alpha_t) * liquid_speed / solid_trans_speed));
            alpha_l = conj(asin(sin(alpha_t) * solid_long_speed / solid_trans_speed));
            
            ct_cl2 = (solid_trans_speed * solid_trans_speed) / (solid_long_speed * solid_long_speed);
            cos_2_alpha_t = cos(2 * alpha_t);
            N = ( ...
                ct_cl2 * sin(2 * alpha_l) .* sin(2 * alpha_t) ...
                + cos_2_alpha_t .* cos_2_alpha_t ...
                + liquid_dens * liquid_speed / (solid_dens * solid_long_speed) * cos(alpha_l) ./ cos(alpha_fluid) ...
            );
            
            if ~mode_out % If output mode is longitudinal.
                % refl_l.
                coeff = -sin(4 * alpha_t) ./ N;
            
                % For conversion to displacement: z = 
                % (rho_inc * vel_inc)/(rho_out * c_out)
                z = (solid_trans_speed) / (solid_long_speed);
                coeff = z * coeff;
                
            else % Output mode is transverse.
                % refl_t.
                coeff = ( ...
                    ct_cl2 * sin(2 * alpha_l) .* sin(2 * alpha_t) ...
                    - cos_2_alpha_t .* cos_2_alpha_t ...
                    - liquid_dens * liquid_speed / (solid_dens * solid_long_speed) * cos(alpha_l) ./ cos(alpha_fluid) ...
                ) / N;
            
                % For conversion to displacement: z = 
                % (rho_inc * vel_inc)/(rho_out * c_out). Note that when reflecting,
                % if inc_mode == out_mode, then z = 1, so no conversion required.
                
            end
        
        end
    else % Incident material must be liquid. Output material will be liquid.
        % Note that this means both input and output mode must be
        % longitudinal.
        alpha_fluid = angle_inc;
        alpha_l = conj(asin(sin(alpha_fluid) * solid_long_speed / liquid_speed));
        alpha_t = conj(asin(sin(alpha_fluid) * solid_trans_speed / liquid_speed));
        
        ct_cl2 = (solid_trans_speed * solid_trans_speed) / (solid_long_speed * solid_long_speed);
        cos_2_alpha_t = cos(2 * alpha_t);
        N = ( ...
            ct_cl2 * sin(2 * alpha_l) .* sin(2 * alpha_t) ...
            + cos_2_alpha_t .* cos_2_alpha_t ...
            + liquid_dens * liquid_speed / (solid_dens * solid_long_speed) * cos(alpha_l) ./ cos(alpha_fluid) ...
        );
    
        % refl.
        coeff = ( ...
            ct_cl2 * sin(2 * alpha_l) .* sin(2 * alpha_t) ...
            + cos_2_alpha_t .* cos_2_alpha_t ...
            - (liquid_dens * liquid_speed * cos(alpha_l)) ./ (solid_dens * solid_long_speed * cos(alpha_fluid)) ...
        ) ./ N;
            
        % For conversion to displacement: z = 
        % (rho_inc * vel_inc)/(rho_out * c_out). Note that when reflecting,
        % if inc_mode == out_mode, then z = 1, so no conversion required.
        
    end
        
end

end
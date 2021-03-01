function scat_info = fn_scat_info(varargin)
% Function comparable to fn_path_info which collects the information passed
% in as input into a single structure for use in calculation.
% INPUTS:
% - type : string
%       String which labels the type of scatterer. Currently allowed values
%       are limited to 'point', 'sdh', 'crack' and 'image'. Type of
%       scatterer will dictate the number of additional input parameters.
% POINT INPUTS:
% - vL : double
%       Longitudinal velocity of the wave in the inspection block where the
%       scatterer is located.
% - vT : double
%       Shear velocity of the wave in the inspection block where the
%       scatterer is located.
% - image_block : array (no_scats, 3)
%       3D coordinates of all of the points in the inspection block. It is not
%       optional as it is required for TFMs, despite the fact that it will
%       be overwritten when computing sensitivity maps. As it is only a
%       small field, this will take very little extra time to write and
%       memory to contain.
% SDH INPUTS:
% - radius : double
%       The radius of the side-drilled hole.
% - lambdaL : array (no_freqs, 1)
%       The wavelengths of the longitudinal wave in the inpsection block
%       where the scatterer is located. If using the single-frequency
%       model, no_freqs = 1. If using the multi-frequency model, the
%       frequencies should be the non-zero values output from the
%       fn_create_input_signal function.
% - lambdaT : array (no_freqs, 1)
%       The wavelengths of the shear wave in the inpsection block where the
%       scatterer is located. If using the single-frequency
%       model, no_freqs = 1. If using the multi-frequency model, the
%       frequencies should be the non-zero values output from the
%       fn_create_input_signal function.
% - angle : double
%       The angle that the side-drilled hole makes in radians with respect
%       to the zero-angle in the solid.
% - image_block : array (no_scats, 3)
%       3D coordinates of all of the side-drilled holds in the inspection
%       block. It is not optional as it is required for TFMs, despite the
%       fact that it will be overwritten when computing sensitivity maps.
%       As it is only a small field, this will take very little extra time 
%       to write and memory to contain.
% - ang_pts_over_2pi : OPTIONAL integer : DEFAULT = 0
%       The number of discrete points for which the scattering amplitude
%       will be precomputed. If default, no precomputation will be
%       performed.
% CRACK INPUTS:
% - vL : double
%       Longitudinal velocity of the wave in the inspection block where the
%       scatterer is located.
% - vT : double
%       Shear velocity of the wave in the inspection block where the
%       scatterer is located.
% - dens : double
%       Density of the inspection material where the scatterer is located.
% - freq : array (no_freqs, 1)
%       The non-zero frequencies to be modelled when calculating the
%       scattering amplitudes. If using the single-frequency
%       model, no_freqs = 1. If using the multi-frequency model, the
%       frequencies should be the non-zero values output from the
%       fn_create_input_signal function.
% - crack_length : double
%       Length of the crack from end to end.
% - angle : double
%       The angle that the crack makes in radians with respect to the zero-
%       angle in the solid.
% - image_block : OPTIONAL array (no_scats, 3)
%       3D coordinates of all of the crack in the inspection block. It is not
%       optional as it is required for TFMs, despite the fact that it will
%       be overwritten when computing sensitivity maps. As it is only a
%       small field, this will take very little extra time to write and
%       memory to contain.
% - nodes_per_wavelength : OPTIONAL integer : DEFAULT = 20
%       The number of nodes per wavelength to model (NEED TO READ INTO WHAT
%       THIS DOES!)
% - ang_pts_over_2pi : OPTIONAL integer : DEFAULT = 0
%       The number of discrete points for which the scattering amplitude
%       will be precomputed. If default, no precomputation will be
%       performed.
% IMAGE INPUTS:
% - image_block : OPTIONAL array (no_scats, 3)
%       3D coordinates of all of the crack in the inspection block. It is not
%       optional as it is required for TFMs, despite the fact that it will
%       be overwritten when computing sensitivity maps. As it is only a
%       small field, this will take very little extra time to write and
%       memory to contain.
%
% OUTPUTS:
% - scat_info : struct (1, 1)
%       Structure containing all of the information of a single scatterer
%       type. Note that multiples of these structures can be provided to
%       the simulation code to simulate multiple types.

if varargin{1} == "point"
    nreqargs = 4;

    % Required arguments
    scat_info.type = varargin{1};
    scat_info.vL = varargin{2};
    scat_info.vT = varargin{3};
    scat_info.image_block = varargin{4};

elseif varargin{1} == "sdh"
    nreqargs = 6;
    
    % Required arguments
    scat_info.type = varargin{1};
    scat_info.radius = varargin{2};
    scat_info.lambdaL = varargin{3};
    scat_info.lambdaT = varargin{4};
    scat_info.angle = varargin{5};
    scat_info.image_block = varargin{6};
    
    % Optional arguments
    if nargin > nreqargs
        i = nreqargs + 1;
        while(i <= size(varargin, 2))
            switch lower(varargin{i})
                case 'ang_pts_over_2pi'
                    ang_pts_over_2pi = varargin{i+1};
                    if ang_pts_over_2pi ~= 0
                        scat_info.matrix = fn_scattering_matrix(scat_info, ang_pts_over_2pi);
                    end
                    i = i+2;
                otherwise
                    error('fn_scat_info: Invalid argument provided for sdh.')
            end
        end
    end

elseif varargin{1} == "crack"
    nreqargs = 8;
    
    % Required arguments
    scat_info.type = varargin{1};
    scat_info.vL = varargin{2};
    scat_info.vT = varargin{3};
    scat_info.dens = varargin{4};
    scat_info.freq = varargin{5};
    scat_info.crack_length = varargin{6};
    scat_info.angle = varargin{7};
    scat_info.image_block = varargin{8};
    
    % Optional arguments
    if nargin > nreqargs
        i = nreqargs + 1;
        while(i <= size(varargin, 2))
            switch lower(varargin{i})
                case 'nodes_per_wavelength'
                    scat_info.nodes_per_wavelength = varargin{i+1};
                    i = i+2;
                case 'ang_pts_over_2pi'
                    ang_pts_over_2pi = varargin{i+1};
                    if ang_pts_over_2pi ~= 0
                        scat_info.matrix = fn_scattering_matrix(scat_info, ang_pts_over_2pi);
                    end
                    i = i+2;
                otherwise
                    error('fn_scat_info: Invalid argument provided for crack.')
            end
        end
    end
    
    if ~isfield(scat_info, 'nodes_per_wavelength')
        scat_info.nodes_per_wavelength = 20;
    end

elseif varargin{1} == "image"
    scat_info.type = varargin{1};
    scat_info.image_block = varargin{2};

else
    error('fn_scat_info: Invalid scatterer type.')
end

end
function scat_info = fn_scat_info(varargin)
% Function comparable to fn_path_info which collects the information passed
% in as input into a single structure for use in calculation.
% INPUTS:
% - type : string
%       String which labels the type of scatterer. Currently allowed values
%       are limited to 'point', 'sdh', 'crack' and 'image'. Type of
%       scatterer will dictate the number of additional input parameters.
% POINT INPUTS:
% - x : array (no_scats, 1)
%       x-coordinates of the points in the inspection block.
% - y : array (no_scats, 1)
%       y-coordinates of the points in the inspection block.
% - z : array (no_scats, 1)
%       z-coordinates of the points in the inspection block.
% - vL : double
%       Longitudinal velocity of the wave in the inspection block where the
%       scatterer is located.
% - vT : double
%       Shear velocity of the wave in the inspection block where the
%       scatterer is located.
% SDH INPUTS:
% - x : array (no_scats, 1)
%       x-coordinates of the side-drilled holes.
% - y : array (no_scats, 1)
%       x-coordinates of the side-drilled holes.
% - z : array (no_scats, 1)
%       x-coordinates of the side-drilled holes.
% - r : array (no_scats, 1)
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
% - ang_pts_over_2pi : OPTIONAL integer : DEFAULT = 0
%       The number of discrete points for which the scattering amplitude
%       will be precomputed. If default, no precomputation will be
%       performed.
% CRACK INPUTS:
% - x : array (no_scats, 1)
%       x-coordinates of the cracks in the inspection block.
% - y : array (no_scats, 1)
%       y-coordinates of the cracks in the inspection block.
% - z : array (no_scats, 1)
%       z-coordinates of the cracks in the inspection block.
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
% - nodes_per_wavelength : OPTIONAL integer : DEFAULT = 20
%       The number of nodes per wavelength to model (NEED TO READ INTO WHAT
%       THIS DOES!)
% - ang_pts_over_2pi : OPTIONAL integer : DEFAULT = 0
%       The number of discrete points for which the scattering amplitude
%       will be precomputed. If default, no precomputation will be
%       performed.
% IMAGE INPUTS:
% - x : array (no_scats, 1)
%       x-coordinates of the points in the inspection block.
% - y : array (no_scats, 1)
%       y-coordinates of the points in the inspection block.
% - z : array (no_scats, 1)
%       z-coordinates of the points in the inspection block.
%
% OUTPUTS:
% - scat_info : struct (1, 1)
%       Structure containing all of the information of a single scatterer
%       type. Note that multiples of these structures can be provided to
%       the simulation code to simulate multiple types.

if varargin{1} == "point"
    nreqargs = 6;

    % Required arguments
    scat_info.type = varargin{1};
    if iscell(varargin{2})
        scat_info.x = cell2mat(varargin{2});
    else
        scat_info.x = varargin{2};
    end
    if iscell(varargin{3})
        scat_info.y = cell2mat(varargin{3});
    else
        scat_info.y = varargin{3};
    end
    if iscell(varargin{4})
        scat_info.z = cell2mat(varargin{4});
    else
        scat_info.z = varargin{4};
    end
    scat_info.vL = varargin{5};
    scat_info.vT = varargin{6};

elseif varargin{1} == "sdh"
    nreqargs = 8;
    
    % Required arguments
    scat_info.type = varargin{1};
    if iscell(varargin{2})
        scat_info.x = cell2mat(varargin{2});
    else
        scat_info.x = varargin{2};
    end
    if iscell(varargin{3})
        scat_info.y = cell2mat(varargin{3});
    else
        scat_info.y = varargin{3};
    end
    if iscell(varargin{4})
        scat_info.z = cell2mat(varargin{4});
    else
        scat_info.z = varargin{4};
    end
    if iscell(varargin{5})
        scat_info.r = cell2mat(varargin{5});
    else
        scat_info.r = varargin{5};
    end
    scat_info.lambdaL = varargin{6};
    scat_info.lambdaT = varargin{7};
    scat_info.angle = varargin{8};
    
    % Optional arguments
    if nargin > nreqargs
        ii = nreqargs + 1;
        while(ii <= size(varargin, 2))
            switch lower(varargin{ii})
                case 'ang_pts_over_2pi'
                    ang_pts_over_2pi = varargin{ii+1};
                    scat_info.ang_pts_over_2pi = ang_pts_over_2pi;
                    if ang_pts_over_2pi ~= 0
                        scat_info.matrix = fn_scattering_matrix(scat_info, ang_pts_over_2pi);
                    end
                    ii = ii+2;
                case 'fmc_mask'
                    fmc_mask = varargin{ii+1};
                    scat_info.fmc_mask = fmc_mask;
                    ii = ii+2;
                otherwise
                    error('fn_scat_info: Invalid argument provided for sdh.')
            end
        end
    else
        scat_info.ang_pts_over_2pi = 0;
    end

elseif varargin{1} == "crack"
    nreqargs = 10;
    
    % Required arguments
    scat_info.type = varargin{1};
    if iscell(varargin{2})
        scat_info.x = cell2mat(varargin{2});
    else
        scat_info.x = varargin{2};
    end
    if iscell(varargin{3})
        scat_info.y = cell2mat(varargin{3});
    else
        scat_info.y = varargin{3};
    end
    if iscell(varargin{4})
        scat_info.z = cell2mat(varargin{4});
    else
        scat_info.z = varargin{4};
    end
    scat_info.vL = varargin{5};
    scat_info.vT = varargin{6};
    scat_info.dens = varargin{7};
    scat_info.freq = varargin{8};
    scat_info.crack_length = varargin{9};
    scat_info.angle = varargin{10};
    
    % Optional arguments
    if nargin > nreqargs
        ii = nreqargs + 1;
        while(ii <= size(varargin, 2))
            switch lower(varargin{ii})
                case 'nodes_per_wavelength'
                    scat_info.nodes_per_wavelength = varargin{ii+1};
                    ii = ii+2;
                case 'ang_pts_over_2pi'
                    ang_pts_over_2pi = varargin{ii+1};
                    if ang_pts_over_2pi ~= 0
                        scat_info.matrix = fn_scattering_matrix(scat_info, ang_pts_over_2pi);
                    end
                    ii = ii+2;
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
    if iscell(varargin{2})
        scat_info.x = cell2mat(varargin{2});
    else
        scat_info.x = varargin{2};
    end
    if iscell(varargin{3})
        scat_info.y = cell2mat(varargin{3});
    else
        scat_info.y = varargin{3};
    end
    if iscell(varargin{4})
        scat_info.z = cell2mat(varargin{4});
    else
        scat_info.z = varargin{4};
    end
    scat_info.fmc_mask = 1;

else
    error('fn_scat_info: Invalid scatterer type.')
end

end
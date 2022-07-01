function Im = fn_create_im(name, xpts, zpts)
% Initialisation function for the image to be calculated. Note that this
% can be used for both images and sensitivity maps.
%
% INPUTS:
% - name : string
%       Simple string which has the name of the view which this image
%       corresponds to.
% - xpts : integer
%       Number of pixels in the x-dimension into which the image can be
%       stored.
% - zpts : integer
%       Number of pixels in the z-dimension into which the image can be
%       stored.
%
% OUTPUTS:
% - Im : struct (1, 1)
%       Initialisation structure with a name and an array to store the
%       image.

Im.name = name;
Im.image = zeros(zpts, xpts);
Im.x = zeros(1, xpts);
Im.z = zeros(zpts, 1);
Im.plotExtras = struct;

end
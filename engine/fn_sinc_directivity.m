function directivity = fn_sinc_directivity(theta, el_width, lambda)
% Computes the directivity from a longitudinal wave. Assumes an immersion
% setup where the waves generated from the array will only be longitudinal
% mode. Note: matlab function defines sinc(x) := sin(pi*x)/(pi*x). Divide
% by pi to account for this, if desired.
%
% INPUTS: 
% - theta : double
%       Outgoing angle from the probe, for which the directivity will be
%       calculated.
% - el_width : double
%       Element width.
% - lambda : double
%       Wavelength of the first leg of the ray.
%
% OUTPUTS : complex
%       Directivity.

% directivity = sinc((el_width/lambda) * sin(theta));
directivity = sinc((el_width/(lambda)) * sin(theta));

end
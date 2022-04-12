function vq = fn_lanczos_interp(x, v, xq, a)
% Lanczos interpolation, similar implementation to arim.
%
% INPUTS:
% - x : array 
%       Sample points
% - v : array
%       Sample values (must have same shape as x)
% - xq : array
%       Query points. May have any shape; vq will return with the same
%       shape as xq.
% - a : int
%       Lanczos interpolation factor
%
% OUTPUTS:
% - vq : array
%       Interpolated values at query points with same shape as xq.

n = length(v);
% xq_idx equiv to lookup_index in arim.
xq_idx = interp1(x, [1:length(v)], xq, 'linear');

i_min = floor(xq_idx) - a + 1;

vq = zeros(size(xq));
for a_idx = 0:2*a - 1
    ii = i_min+a_idx;
    vq = vq + ...
        v(mod(ii, n)) .* sinc(xq_idx - ii) .* sinc((xq_idx - ii)/a);
end

end
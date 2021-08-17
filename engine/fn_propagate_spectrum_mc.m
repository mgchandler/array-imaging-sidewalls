function out_freq_spec = fn_propagate_spectrum_mc(freq, in_freq_spec, tof, varargin);
%USAGE
%	out_freq_spec = fn_propagate_spectrum(freq, in_freq_spec, ph_vel, dists [, amps, use_gpu_if_present]);
%AUTHOR
%	Paul Wilcox (Oct 2007, updated Mar 2012 with GPU support and better CPU version. Old version still available as fn_propagate_spectrum_old)
%SUMMARY
%	Applies phase shifts to spectrum (in_spec) to simulate propagation delays associated
%with one or more propagation distances and returns a phase-shifted spectrum for each
%propagation distance specified.
%	Note that no checks are made to see if propagation distance(s) is/are small enough to be
%accommodated by the number of points in in_freq_spec. If too large a propagation distance is
%specified then signals will be wrapped when out_freq_spec are converted back to time domain.
%INPUTS
%	freq - frequency vector corresponding to spectra
%	in_freq_spec - the spectrum that will be propagated
%	ph_vel - phase velocity: either a single number for non-dispersive propagation 
%				or a 2 column matrix where 1st column is frequency and second is phase velocity
%	dists - required propagation distance or vector if spectra for multiple distances are required
%	amps[ones(size(dists))] - optional vector of corresponding amplitudes
%   use_gpu_if_present[1] - optional binary value
%OUTPUTS
%	out_freq_spec - vector (or matrix if length(dists)>1) of the resulting
%   propagated spectra. If in matrix form, the spectra are in columns, and the number of columns is
%   equal to length(dists). Note that if using the GPU version, this is
%   returned as a gpuArray. You can use gather to convert it back.
%NOTES
%   If using GPU version, keep the output as gpuArray if you are calling this
%   function multiple times (e.g. in a loop) and accumulating the result:
%       for ii = 1:1000
%           out_freq_spec = out_freq_spec + fn_propagate_spectrum(...)
%       end
%   as this is more efficient. Only use gather at the end outside the loop:
%       out_freq_spec = gather(out_freq_spec)
%--------------------------------------------------------------------------

if nargin > 4
   amps = varargin{1};
else
   amps = ones(size(dists));
end;

if nargin > 5
    use_gpu_if_present = varargin{2};
else
    use_gpu_if_present = 1;
end

%get wavenumber spectrum (interpolate if nesc)
omega = 2 * pi * freq;
%do the propagation

%OLD METHOD WITH LOOP
% out_freq_spec = zeros(length(freq),length(dists));
% for ii = 1:length(dists)
%    out_freq_spec(:,ii) = amps(ii) * in_freq_spec(:) .* exp(-i * k(:) * dists(ii));
% end;

%METHOD WITHOUT LOOP (about 30% faster than old method) 
% out_freq_spec = (in_freq_spec(:) * amps(:).') .* exp(-i * k(:) * dists(:).');

%ANOTHER METHOD WITHOUT LOOP (about 34% faster than old method) 
% out_freq_spec = spdiags(in_freq_spec(:), 0, m, m) * exp(-i * k(:) * dists(:).') * spdiags(amps(:), 0, n, n);

%GPU METHOD
if use_gpu_if_present && (exist('gpuDeviceCount') == 2) && (gpuDeviceCount > 0)
    out_freq_spec = (gpuArray(in_freq_spec(:)) * gpuArray(amps(:).')) .* exp(-1i * gpuArray(omega(:)) * gpuArray(tof(:).'));
else
    m = length(omega);
    n = length(amps);
    diag1 = spdiags(in_freq_spec(:), 0, m, m);
    delta = exp(-1i * omega(:) * tof(:).');
    diag2 = spdiags(amps(:), 0, n, n);
    out_freq_spec = spdiags(in_freq_spec(:), 0, m, m) * exp(-1i * omega(:) * tof(:).') * spdiags(amps(:), 0, n, n);
end
return;
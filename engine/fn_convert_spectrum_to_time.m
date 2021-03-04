function [time, out_time_sig] = fn_convert_spectrum_to_time(freq, out_freq_spec, fft_pts, time_pts, varargin);
%USAGE
%	[time, out_time_sig] = fn_convert_spectrum_to_time(freq, out_freq_spec, fft_pts, time_pts [,return_hilbert]);
%AUTHOR
%	Paul Wilcox (2007)
%SUMMARY
%	Converts one or more complex frequency spectra to real time domain signals.
%Use in conjunction with create_input_signal and propagate_spectrum functions for
%simulating signals.
%INPUTS
%	freq - vector of frequencies, should be from zero to the Nyquist frequency.
%	out_freq_spec - vector or matrix of spectrum or spectra to be transformed. Length
%of vector or number of rows in matrix should be equal to the length(freq).
%	fft_pts - number of points used for Fourier transform when spectrum
%or spectra were initially created.
%	time_pts - number of points required in time domain signals. Must be less than
%or equal to fft_pts.
%   return_hilbert[0] - returns complex hilbert transform of signal is 1.
%OUTPUTS
%	time - vector of times.
%	out_time_sig - vector or matrix of time domain signals
%NOTES
%	Designed for use in conjunction with
%		fn_create_input_signal and fn_propagate_spectrum

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<5
    return_hilbert = 0;
else
    return_hilbert = varargin{1};
end;

if time_pts > fft_pts
   time_pts = fft_pts;
end;

if return_hilbert
    out_time_sig = ifft(out_freq_spec,fft_pts) * 2;
else
    out_time_sig = real(ifft(out_freq_spec,fft_pts)) * 2;
end;

out_time_sig = out_time_sig(1:time_pts,:);

time_step = 1 / (fft_pts*abs(freq(2)-freq(1)));
time = [0:time_pts-1]*time_step;

return;

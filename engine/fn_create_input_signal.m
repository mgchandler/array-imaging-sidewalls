function [time, in_time_sig, freq, in_freq_spec, fft_pts] = fn_create_input_signal(time_pts, centre_freq, time_step, varargin);
%USAGE
%	[time, in_time_sig, freq, in_freq_spec, fft_pts] = fn_create_input_signal(time_pts, centre_freq, time_step [, no_cycles, window_type, centre_time]);
%AUTHOR
%	Paul Wilcox (Oct 2007)
%SUMMARY
%	Creates input time signal and spectrum. Designed for general simulation
%	of standard NDT input signals but also for wave propagation simulation
%	in conjunction with fn_propagate_spectrum and
%	fn_convert_spectrum_to_time.
%   Note that there is no checking of parameters' compatibility.
%INPUTS
%	time_pts - number of points in time signal
%	centre_freq - centre frequency
%	time_step - time step
%	no_cycles [0] - number of cycles (anything < 1 gives pulse)
%	window_type ['gaussian 40'] - window type and necessary parameter string (default = 'gaussian 40')
%				'gaussian n' - gaussian based on n dB down points
%				'hanning' - hanning
%	centre_time [0.0] - time of centre of signal
%OUTPUTS
%	time - real vector of times, length equal to time_pts
%	in_time_sig - real vector of the generated sime signal, length equal to time_pts. 
%	freq - real vector of frequencies starting at zero and finishing at the Nyquist frequency. Length is
%equal to fft_pts/2+1
%	in_freq_spec - complex vector of the positive half of the frequency spectrum of time_sig. Length is
%equal to fft_pts/2+1
%	fft_pts - the number of points used in the Fourier transform to calculate the freq_spec
%NOTES
%	This function is specifically designed for use in conjunction with
%		out_freq_spec=propagate_spectrum(freq, in_freq_spec, ph_vel, dists, amps);
%		[time, out_time_sig] = convert_spectrum_to_time(freq, out_freq_spec, fft_pts, time_pts);
%	for simulating time-domain signals
%	Note that the output in_time_sig is not used for the simulation of time-trace with propagate_spectrum.
%	Its main use is for checking that this function works. 
%	Also note that if centre_time is zero (the default) then only the second half of the signal will be
%	visible in in_time_sig. Don't worry, all the information is in freq_spec.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

default_no_cycles = 0;
default_window_type = 'gaussian 40';
default_centre_time = 0;

if nargin>3
   no_cycles = varargin{1};
else
   no_cycles = default_no_cycles;
end;

if nargin>4
   window_type = varargin{2};
else
   window_type = default_window_type;
end;

if nargin>5
   centre_time = varargin{3};
else
   centre_time = default_centre_time;
end;

%%%do it
time = (([1:time_pts]-1)*time_step)';
tmid = max(time) / 2;
tmax = max(time);

if centre_time < 0
   tx = 0;
end;

if centre_time > tmax
   tx = tmax;
end;

if no_cycles < 1.0
   carrier = ones(time_pts, 1);
   no_cycles = 0.5;
else
	carrier = sin(2*pi*centre_freq*(time - tmid));
end;

if strncmp(window_type,'gaussian',8)
   db = sscanf(window_type,'gaussian %f');
   window = fn_gaussian(time_pts, 0.5, no_cycles / centre_freq / tmax / 2,db);
   time_sig = carrier .* window;
end;

if strncmp(window_type,'hanning',7)
   window = fn_hanning(time_pts, 0.5, no_cycles / centre_freq / tmax / 2);
   time_sig = carrier .* window;
end;

%do spectrum
duration = no_cycles / centre_freq;
min_pts = ceil((tmax + duration / 2)/time_step);
fft_pts = 2 ^ nextpow2(min_pts);
fstep = 1/(fft_pts * time_step);
freq = ((0:fft_pts / 2) * fstep)';
in_freq_spec = fft(time_sig, fft_pts);
in_freq_spec = in_freq_spec(1:fft_pts / 2 + 1);
in_freq_spec = in_freq_spec .* exp(2 * pi * i * freq * (tmid - centre_time));
in_time_sig = real(ifft(in_freq_spec, fft_pts)) * 2;
in_time_sig = in_time_sig(1:time_pts);
sf = 1 / max(abs(in_time_sig));
in_time_sig = in_time_sig * sf;
in_freq_spec = in_freq_spec * sf;

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function window = fn_hanning(number_of_points, peak_pos_fract, half_width_fract, varargin);

x = linspace(0, 1, number_of_points)';
y = 0.5 * (1 + cos((x - peak_pos_fract) / half_width_fract * pi));
window = y .* ((x >= (peak_pos_fract - half_width_fract)) & (x <= (peak_pos_fract + half_width_fract)));

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function window=fn_gaussian(number_of_points, peak_pos_fract, half_width_fract, varargin);

default_db_down = 40;
default_force_zero = 0;
if nargin < 4;
   db_down = default_db_down;
else
   db_down = varargin{1};
end;
if nargin < 5;
   force_zero = default_force_zero;
else
   force_zero = varargin{2};
end;

fract = 10.0 ^ (-db_down / 20.0);
r = (linspace(0, 1, number_of_points) - peak_pos_fract)';
r1 = half_width_fract / ((-log(fract)) ^ 0.5);
window = exp(-(r / r1) .^ 2);
if force_zero
	window(find(window < fract)) = 0.0;
end;

return;
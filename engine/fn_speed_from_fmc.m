function v_L = fn_speed_from_fmc(FMC_time, FMC_time_data, tx, rx, backwall_offset, varargin)
% Measures the wave velocity of the fastest wave speed using reverberations
% from a back wall. Assumes that the array is in contact with the sample,
% and that geometrical artefacts are present (i.e. a non-defective FMC has
% not been subtracted for FE).
%
% INPUTS:
%   - FMC_time : double array
%       Shape (Nt, 1). Time domain points.
%   - FMC_time_data : double/complex array
%       Shape (Nt, Nprobes^2). Actual FMC data. Can be either double or
%       complex as abs() is called.
%   - tx and rx : double array
%       Shape (1, Nprobes^2). Tx and Rx indices of each term in
%       FMC_time_data.
%   - backwall_offset : double
%       Distance to the back wall from the array.
%   - frac_of_max : double
%       Fraction of the largest peak acting as a threshold, above which
%       signals will be considered as reverberations.
%
% NOTES:
%   This function was designed to extract longitudinal velocity from FE
%   FMC data modeling an aluminium sample, for which vS ~= .5 * vL, thus it
%   is assumed that reverberation signals are quite clean. If there are
%   several highly prominent peaks clustered together, this will likely
%   give poorer results.

assert(size(FMC_time, 1) == size(FMC_time_data, 1), "fn_speed_from_fmc: FMC_time and FMC_time_data have incompatible number of time points.")
assert(size(FMC_time_data, 2) == size(tx, 2), "fn_speed_from_fmc: FMC_time_data and tx have incompatible number of timetraces.")
assert(all(size(tx) == size(rx)), "fn_speed_from_fmc: tx and rx are not the same size.")
if size(FMC_time, 2) ~= 1
    FMC_time = FMC_time(:, 1);
end

if nargin > 5
    frac_of_max = varargin{1};
else
    frac_of_max = 20;
end

truncated_FMC = sum(abs(FMC_time_data(:, tx==rx)), 2);
[~, locs, ~, p] = findpeaks(truncated_FMC);

times = sort(FMC_time(locs(p>max(p)/frac_of_max)));
tof = mean(diff(times));
std_tof = std(diff(times));

if tof/std_tof < 200
    warning("fn_speed_from_fmc: Mean-to-StDev ratio for reverberation time is too small (avg/std = %.3f). FMC likely too noisy.", tof/std_tof)
end

v_L = 2*backwall_offset / tof;

end
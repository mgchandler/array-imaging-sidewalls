function fn_plot_FMC_at_time(FMC_data, FMC_time, path_info_tx, path_info_rx, pixel_coords, savename)
% Plots a grayscale FMC. The ToF to a specific pixel in a TFM is
% calculated, and the part of the FMC which contributes to this pixel for
% each tx-rx pair is plotted in colour on top of the grayscale image.
%
% Note that for each tx-rx pair, a band of times are plotted in colour
% which is centred on the ToF. For very fine time sampling, plotting a
% single time point is meaningless thus a band makes the resulting data
% easier to interpret.
%
% INPUTS:
% - FMC_data : complex OR double array (num_time_pts, num_els^2)
%       FMC data obtained either from simulation or experimental data. If
%       type is complex, it will be converted to double using the abs()
%       method.
% - FMC_time : double array (num_time_pts, 1)
%       Corresponding time scale.
% - path_info_tx : struct (1, 1)
%       Output of the fn_path_info function. Used as the transmit path of
%       the view when calculating ToF to the pixel of choice.
% - path_info_rx : struct (1, 1)
%       Output of the fn_path_info function. Used as the transmit path of
%       the view when calculating ToF to the pixel of choice.
% - pixel_coord : double array (1, 3)
%       Real coordinate of the pixel for which the ToF will be calculated.
% - savename : string
%       Name which the figure generated will be saved as. If it's an empty
%       string, no figure will be saved.

probe_els_2 = size(FMC_data, 2);
time_pts  = size(FMC_time, 1);
num_scats = size(pixel_coords, 1);

% Fraction of the total time domain which should be plotted in colour for
% each pixel. Note that if multiple pixels in the TFM are desired, this
% value should be smaller.
perc = .002;

assert(size(FMC_data, 1) == time_pts, "fn_plot_FMC_at_time: FMC_data and FMC_time require same dimension in axis=1.")
% assert(num_scats == 1, "fn_plot_FMC_at_time: Incorrect number of pixels provided.")
assert(size(pixel_coords, 2) == 3, "fn_plot_FMC_at_time: Incorrect number of dimensions in pixel.")

% Work towards getting the min_times to each pixel from each element
pixel_info = fn_scat_info("image", pixel_coords(:, 1), pixel_coords(:, 2), pixel_coords(:, 3));
tx_path = fn_compute_ray(pixel_info, path_info_tx);
rx_path = fn_compute_ray(pixel_info, path_info_rx);
view = fn_create_view(tx_path, rx_path);
% Collapse min_times into one vector. view.min_times has shape
% (probe_els^2, num_scats) so vector element (1) should have probe=1
% scat=1; element (2): probe=2 scat=1; element (probe_els+1): probe=1
% scat=2 etc.
min_times = view.min_times(:);

% Find where in FMC_time the ToF to each pixel is. Restore the original
% shape of these indices.
[row, col] = find(abs(min_times - FMC_time.') == min(abs(min_times - FMC_time.'), [], 2));
[~,I] = sort(row);
idxs = col(I);
idxs = reshape(idxs, size(view.min_times));

% Make a grayscale FMC
lower_limit = 1/4;
gray_data = abs(FMC_data')/max(abs(FMC_data(:))) * (1-lower_limit) + lower_limit;
grayFMC = repmat(gray_data, [1, 1, 3]);

FMC_data_t = FMC_data.';

% Assign colour to the bits of interest.
color_FMC = zeros(size(FMC_data.'));
for scat = 1:num_scats
    color_FMC(abs(view.min_times(:, scat) - FMC_time.') < perc*FMC_time(end)) = ...
        abs(FMC_data_t(abs(view.min_times(:, scat) - FMC_time.') < perc*FMC_time(end))) ./ ...
        max(abs(FMC_data(:)));
end

figure(2);
t = tiledlayout(1,1);%(1,1,'Padding','tight');
% t.Units = 'inches';
% t.OuterPosition = [0.15 0.15 15 9];
nexttile;
imagesc(FMC_time*10^6, 1:probe_els_2, grayFMC);
hold on
fg = imagesc(FMC_time*10^6, 1:probe_els_2, color_FMC);
set(fg, 'AlphaData', im2double(color_FMC~=0))
title(sprintf("view %s", view.name))
xlabel('Time (us)')
ylabel('tx-rx Index')
if savename ~= ""
    if contains(savename, ".fig")
        savefig(savename)
    else
        exportgraphics(t, savename, 'Resolution', 500)
    end
end

end
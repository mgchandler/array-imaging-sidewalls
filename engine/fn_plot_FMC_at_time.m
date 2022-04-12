function fn_plot_FMC_at_time(FMC_data, FMC_time, path_info_tx, path_info_rx, pixel_coord, savename)
% Plots a grayscale FMC. The ToF to a specific pixel in a TFM is
% calculated, and the part of the FMC which contributes to this pixel for
% each tx-rx pair is plotted in colour on top of the grayscale image.
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

assert(size(FMC_data, 1) == size(FMC_time, 1), "fn_plot_FMC_at_time: FMC_data and FMC_time require same dimension in axis=1.")
assert(size(pixel_coord, 1) == 1, "fn_plot_FMC_at_time: Incorrect number of pixels provided.")
assert(size(pixel_coord, 2) == 3, "fn_plot_FMC_at_time: Incorrect number of dimensions in pixel.")

pixel_info = fn_scat_info("image", pixel_coord(1), pixel_coord(2), pixel_coord(3));
tx_path = fn_compute_ray(pixel_info, path_info_tx);
rx_path = fn_compute_ray(pixel_info, path_info_rx);

view = fn_create_view(tx_path, rx_path);

[row, col] = find(abs(view.min_times - FMC_time') == min(abs(view.min_times - FMC_time'), [], 2));
[~,I] = sort(row);
idxs = col(I);

time_idxs = zeros(size(FMC_data, 2), 2);
time_idxs(:, 1) = idxs - round(0.005 * length(FMC_time));
time_idxs(:, 2) = idxs + round(0.005 * length(FMC_time));

lower_limit = 1/4;

gray_data = abs(FMC_data')/max(abs(FMC_data(:))) * (1-lower_limit) + lower_limit;
grayFMC = repmat(gray_data, [1, 1, 3]);

color_FMC = abs(FMC_data')/max(abs(FMC_data(:)));
for ii = 1:size(time_idxs, 1)
    color_FMC(ii, 1:time_idxs(ii, 1)) = 0;
    color_FMC(ii, time_idxs(ii, 2):end) = 0;
end

figure(2);
t = tiledlayout(1,1,'Padding','tight');
t.Units = 'inches';
t.OuterPosition = [0.15 0.15 5 3];
nexttile;
imagesc(FMC_time*10^6, [1:size(FMC_data, 2)], grayFMC);
hold on
fg = imagesc(FMC_time*10^6, [1:size(FMC_data, 2)], color_FMC);
set(fg, 'AlphaData', im2double(color_FMC~=0))
title(sprintf("%s, view %s highlighted", savename(1:3), view.name))
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
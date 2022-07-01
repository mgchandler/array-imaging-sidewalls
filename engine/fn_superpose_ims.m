function fn_superpose_ims(Im1, Im2, varargin)
% Produces subplots which superpose Im2 onto Im1. It is assumed that Im2
% will contain a higher resolution subset of Im1, as Im2 will be plotted on
% top of Im1.
%
% INPUTS:
%   - x1, z1, x2, z2 : double array
%       x- and z- axes for Im1 and Im2 respectively.
%   - Im1, Im2 : struct array
%       Expected to have form output from fn_create_im(). "image" fields
%       are expected to be on the same scale.
%   - morePlots : struct array
%       Contains additional features to plot on top of the image. Intended
%       to allow geometry, array elements and ray paths to be plotted as
%       well if desired.

%% Argument validity checks

UC = 1e3;
db_range_for_output = 40;

number_of_ims = size(Im1, 1);
if number_of_ims == 3
    plot_x = 3;
    plot_z = 1;
elseif number_of_ims == 21
    plot_x = 3;
    plot_z = 7;
elseif number_of_ims == 55
    plot_x = 11;
    plot_z = 5;
else
    assert(any(number_of_ims == [3, 21, 55]), "fn_superpose_ims: incorrect number of images. Expected 3 (direct), 21 (full-skip) or 55 (double-skip).")
end
assert(number_of_ims == size(Im2, 1), "fn_superpose_ims: incompatible number of images")
for im = 1:number_of_ims
    assert(strcmp(Im1(im).name, Im2(im).name), sprintf("fn_superpose_ims: image %d are not the same.", im))
end
if nargin > 6
    morePlots = varargin{1};
    assert(isstruct(morePlots), "fn_superpose_ims: morePlots expected to be a struct.")
else
    morePlots = 0;
end

%% Setup

% Get max signal across all images.
max_ = 0;
for im = 1 : number_of_ims
    if max(abs(Im1(im).image(:))) > max_
        max_ = max(abs(Im1(im).image(:)));
    end
    if max(abs(Im2(im).image(:))) > max_
        max_ = max(abs(Im2(im).image(:)));
    end
end

% Re-compute db_image based on new max_.
for im = 1 : number_of_ims
   Im1(im).db_image = 20 * log10(abs(Im1(im).image) ./ max_); 
   Im2(im).db_image = 20 * log10(abs(Im2(im).image) ./ max_); 
end

%% Plotting

fig = figure(1);
t = tiledlayout(plot_z, plot_x, 'TileSpacing', 'Compact');
for im = 1:number_of_ims
%     im = im_idxs(im1);
    h(im) = nexttile;
    imagesc(Im1(im).x*UC, Im1(im).z*UC, Im1(im).db_image);
    hold on
    imagesc(Im2(im).x*UC, Im2(im).z*UC, Im2(im).db_image);
    title(Im1(im).name)
    caxis([-db_range_for_output, 0])
    
    for item = 1:size(Im1(im).plotExtras, 2)
        plot(Im1(im).plotExtras(item).x*UC, Im1(im).plotExtras(item).z*UC, 'Color', Im1(im).plotExtras(item).color, 'LineStyle', Im1(im).plotExtras(item).lineStyle, 'Marker', Im1(im).plotExtras(item).marker)
    end
    
%     plot(probe_coords(:, 1)*UC, probe_coords(:, 3)*UC, 'go');
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for wall = 1:size(geometry, 1)
%         plot(geometry(wall).coords(:, 1)*UC, geometry(wall).coords(:, 3)*UC, 'r')
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if boxsize ~= 0
%         for s = 1 : size(scat_info.x, 1)
%             if ~strcmp(scat_info.type, 'image')
%                 new_box_x = scat_info.x(s) + scat_info.r(s)/2*(sin(mean(Views(im).scat_inc_angles))+sin(mean(Views(im).scat_out_angles)));
%                 new_box_z = scat_info.z(s) + scat_info.r(s)/2*(cos(mean(Views(im).scat_inc_angles))+cos(mean(Views(im).scat_out_angles)));
%                 rectangle('Position', [scat_info.x(s)*UC - boxsize*UC/2, scat_info.z(s)*UC - boxsize*UC/2, boxsize*UC, boxsize*UC], 'EdgeColor', 'r');
%                 rectangle('Position', [new_box_x*UC - boxsize*UC/2, new_box_z*UC - boxsize*UC/2, boxsize*UC, boxsize*UC], 'EdgeColor', 'g');
% 
%                 for leg = 1:size(Views(im).path_1.coords, 3)-1
%                     plot([mean(Views(im).path_1.coords(:, s, leg, 1))*UC, mean(Views(im).path_1.coords(:, s, leg+1, 1))*UC], ...
%                          [mean(Views(im).path_1.coords(:, s, leg, 3))*UC, mean(Views(im).path_1.coords(:, s, leg+1, 3))*UC], ...
%                     'Color', [.5,.5,.5])
%                 end
%                 for leg = 1:size(Views(im).path_2.coords, 3)-1
%                     plot([mean(Views(im).path_2.coords(:, s, leg, 1))*UC, mean(Views(im).path_2.coords(:, s, leg+1, 1))*UC], ...
%                          [mean(Views(im).path_2.coords(:, s, leg, 3))*UC, mean(Views(im).path_2.coords(:, s, leg+1, 3))*UC], ...
%                     'Color', [.5,.5,.5])
%                 end
%             end
%         end
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if mod(im, plot_x) ~= 1
        set(gca, 'yticklabel', {[]})
    end
    if im <= number_of_ims - plot_x
        set(gca, 'xticklabel', {[]})
    end
    
    axis equal; axis tight;
    xlim([min(min(Im1(im).x), min(Im2(im).x))*UC, max(max(Im1(im).x), max(Im2(im).x))*UC])
    ylim([min(min(Im1(im).z), min(Im2(im).z))*UC, max(max(Im1(im).z), max(Im2(im).z))*UC])
end
xlabel(t, 'x (mm)')
ylabel(t, 'z (mm)')

c = colorbar(h(1), 'AxisLocation','in');
c.Layout.Tile = 'north';
c.Label.String = 'dB';

end
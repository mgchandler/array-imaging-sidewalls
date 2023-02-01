function fn_image_from_mat(Ims)
% Plot all the Ims, expected to be output from fn_image_tfm().

UC = 1e3;
db_range_for_output = 40;
Number_of_ims = length(Ims);
if Number_of_ims == 1
    plot_x = 1;
    plot_z = 1;
elseif Number_of_ims == 3
    plot_x = 3;
    plot_z = 1;
elseif Number_of_ims == 21
    plot_x = 3;
    plot_z = 7;
elseif Number_of_ims == 55
    plot_x = 11;
    plot_z = 5;
else
    error('fn_sens: Unexpected number of images being plotted.\n%d image(s) being plotted.', Number_of_ims)
end

fig = figure;
t = tiledlayout(plot_z, plot_x, 'TileSpacing', 'Compact');
for im = 1:Number_of_ims
%     im = im_idxs(im1);
    h(im) = nexttile;
    g = imagesc(Ims(im).x*UC, Ims(im).z*UC, Ims(im).db_image);
    hold on
    title(Ims(im).name)
    caxis([-db_range_for_output, 0])
    if isfield(Ims(im), 'plotExtras')
        for ex = 1:length(Ims(im).plotExtras)
            plot(Ims(im).plotExtras(ex).x*UC, Ims(im).plotExtras(ex).z*UC, ...
                 'Color', Ims(im).plotExtras(ex).color, ...
                 'LineStyle', Ims(im).plotExtras(ex).lineStyle, ...
                 'Marker', Ims(im).plotExtras(ex).marker ...
            )
        end
    end
    
    if and(mod(im+1, plot_x)~=0, Number_of_ims~=1)
        set(gca, 'yticklabel', {[]})
    end
    if im <= Number_of_ims - plot_x
        set(gca, 'xticklabel', {[]})
    end
    
    set(g, 'AlphaData', ~or(isnan(Ims(im).db_image), isinf(Ims(im).db_image)))
    
    axis equal; axis tight;
%     xlim([min(Ims(im).x(:))*UC, max(Ims(im).x(:))*UC])
%     ylim([min(Ims(im).z(:))*UC, max(Ims(im).z(:))*UC])
end
% xlabel(t, '$x$ (mm)');
% ylabel(t, '$z$ (mm)');

c = colorbar(h(1), 'AxisLocation','in');
c.Layout.Tile = 'north';
c.Label.String = 'dB';

set(findall(gcf, '-property', 'Fontname'), 'Fontname', 'Serif')
% t.XLabel.Interpreter = 'latex';
% t.YLabel.Interpreter = 'latex';

end
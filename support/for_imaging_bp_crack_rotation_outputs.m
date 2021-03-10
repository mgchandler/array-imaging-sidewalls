cd("C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\array-imaging-sidewalls matlab\bp-output\Crack Rotation - Contact Sensitivity\figs\new")
dd = dir('*.fig');
filenames = {dd.name};

xwidth = 100.0e-3;
N = 100;
half_probe_width = 1.00e-3 * 32 / 2;
xmins = linspace(-xwidth+half_probe_width, -half_probe_width, N);
xmaxs = linspace(+half_probe_width, xwidth-half_probe_width, N);

openname = filenames{1};
savename = openname(1:end-5);

Number_of_ims = 55;

circr = @(radius,rad_ang)  [radius*cos(rad_ang);  radius*sin(rad_ang)];

for f = 1:1%size(filenames, 2)
    openname = filenames{f};
    savename = sprintf('Reformatted %s.png', openname(1:end-4));
    title = sprintf('Crack Rotation θ = %5.1f°', rad2deg(str2double(openname(28:end-4))));
    s(f).name = savename;
    openfig(openname);
    
    fig = gcf;
    h=findobj(gcf,'type','axes');
    figlist = fig.Children;
    data = figlist.Children;
    c=findobj(gcf,'type','colorbar');
    cf = get(c, 'Position');
    set(c,'Position',[cf(1), cf(2)-60, cf(3), cf(4)])
    
    newfig = figure;
    tcl = tiledlayout(11,5);
    for i = 1:numel(figlist)
        figure(figlist(i))
        ax = gca;
        ax.Parent=tcl;
        ax.Layout.Tile=i;
    end
    
    for k=1:Number_of_ims
    	f=get(h(Number_of_ims - k + 1),'Position');
        set(h(Number_of_ims - k + 1),'Position',[f(1), f(2)-30, f(3), f(4)])
%         fprintf('%d: [%.1f, %.1f, %.1f, %.1f]\n', Number_of_ims - k + 1, f(1), f(2), f(3), f(4))
    end
    
    ax = h(1);
    axis manual
    Crack_coords = zeros(2,2);
    Crack_coords(:, 1) = [0, 0];
    Crack_coords(:, 2) = [-15, 15];
    Norm_coords = zeros(2, 2);
    Norm_coords(:, 1) = [0, 15];
    Norm_coords(:, 2) = [0, 0];
    
    
    
    x_shift = 260;
    z_shift = 1050;
    
    ang = str2double(openname(28:end-4));
    Rot_mat = [[cos(ang), -sin(ang)]; [sin(ang), cos(ang)];];
    Crack_coords = Crack_coords * Rot_mat;
    Norm_coords = Norm_coords * Rot_mat;
    
    Arc_coords = circr(10, linspace(pi/2, pi/2-ang, 10)).';
    Arc_avg_x = mean(Arc_coords(:, 1))*2;
    Arc_avg_z = mean(Arc_coords(:, 2))*2;
    
    Crack_coords(:, 1) = Crack_coords(:, 1) - x_shift;
    Crack_coords(:, 2) = Crack_coords(:, 2) - z_shift;
    Norm_coords(:, 1) = Norm_coords(:, 1) - x_shift;
    Norm_coords(:, 2) = Norm_coords(:, 2) - z_shift;
    Arc_coords(:, 1) = Arc_coords(:, 1) - x_shift;
    Arc_coords(:, 2) = Arc_coords(:, 2) - z_shift;
    Arc_avg_x = Arc_avg_x - x_shift;
    Arc_avg_z = Arc_avg_z - z_shift;
    
    hLine_pn = line([-x_shift, -x_shift, -x_shift+1, -x_shift, -x_shift-2], [-z_shift, -z_shift+15, -z_shift+12, -z_shift+15, -z_shift+12], ...
                 'Color', 'k', ...
                 'Clipping', 'off', ...
                 'Linewidth', 1);
    
    hLine_arc = line(Arc_coords(:, 1), Arc_coords(:, 2), ...
                 'Color', 'g', ...
                 'Clipping', 'off', ...
                 'Linewidth', 2);
    
    hLine_crack = line(Crack_coords(:, 1), Crack_coords(:, 2), ...
                 'Color', 'r', ...
                 'Clipping', 'off', ...
                 'Linewidth', 3);
    
    th = text(Arc_avg_x, Arc_avg_z, 'θ');
    
    sgtitle(title)
%     saveas(gcf, savename);
%     close(gcf);
end
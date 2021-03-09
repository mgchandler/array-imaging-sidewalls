cd("C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\array-imaging-sidewalls matlab\bp-output\SDH Scan - Contact Sensitivity\figs")
dd = dir('*.fig');
filenames = {dd.name};

xwidth = 100.0e-3;
N = 100;
half_probe_width = 1.00e-3 * 16;
xmins = linspace(-half_probe_width, -xwidth, N);
xmaxs = linspace(xwidth, half_probe_width, N);

openname = filenames{1};
savename = openname(1:end-5);

Number_of_ims = 21;

for f = 1:size(filenames, 2)
    openname = filenames{f};
    savename = sprintf('%s.png', openname(1:end-4));
    sw_loc = xmaxs(f)-half_probe_width;
    title = sprintf('SDH Sens - Sidewall Position = %4.1f mm', sw_loc*1e3);
    s(f).name = savename;
    openfig(openname);
    
    fig = gcf;
    h=findobj(gcf,'type','axes');
    c=findobj(gcf,'type','colorbar');
    cf = get(c, 'Position');
    set(c,'Position',[cf(1), cf(2)-60, cf(3), cf(4)])
    for k=1:Number_of_ims
    	f=get(h(Number_of_ims - k + 1),'Position');
        set(h(Number_of_ims - k + 1),'Position',[f(1), f(2)-60, f(3), f(4)])
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
    
    
    
%     ang = str2double(openname(28:end-4));
%     Rot_mat = [[cos(ang), -sin(ang)]; [sin(ang), cos(ang)];];
%     Crack_coords = Crack_coords * Rot_mat;
%     Norm_coords = Norm_coords * Rot_mat;
%     
%     Arc_coords = circr(10, linspace(pi/2, pi/2-ang, 10)).';
%     Arc_avg_x = mean(Arc_coords(:, 1))*2;
%     Arc_avg_z = mean(Arc_coords(:, 2))*2;
%     
%     Crack_coords(:, 1) = Crack_coords(:, 1) - 110;
%     Crack_coords(:, 2) = Crack_coords(:, 2) - 610;
%     Norm_coords(:, 1) = Norm_coords(:, 1) - 110;
%     Norm_coords(:, 2) = Norm_coords(:, 2) - 610;
%     Arc_coords(:, 1) = Arc_coords(:, 1) - 110;
%     Arc_coords(:, 2) = Arc_coords(:, 2) - 610;
%     Arc_avg_x = Arc_avg_x - 110;
%     Arc_avg_z = Arc_avg_z - 610;
%     
%     hLine_pn = line([-110, -110, -108, -110, -112], [-610, -595, -598, -595, -598], ...
%                  'Color', 'k', ...
%                  'Clipping', 'off', ...
%                  'Linewidth', 1);
%     
%     hLine_arc = line(Arc_coords(:, 1), Arc_coords(:, 2), ...
%                  'Color', 'g', ...
%                  'Clipping', 'off', ...
%                  'Linewidth', 2);
%     
%     hLine_crack = line(Crack_coords(:, 1), Crack_coords(:, 2), ...
%                  'Color', 'r', ...
%                  'Clipping', 'off', ...
%                  'Linewidth', 3);
%     
%     th = text(Arc_avg_x, Arc_avg_z, 'θ');
    
    sgtitle(title)
    saveas(gcf, savename);
    close(gcf);
end
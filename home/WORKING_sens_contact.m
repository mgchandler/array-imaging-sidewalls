cd('C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\1 - sensitivity\Figures\Contact\Sensitivity')

PITCH = [2.00e-3, 1.00e-3, 0.75e-3, 0.50e-3];
PIXEL = [1.00e-3, 0.50e-3];
WALLS = [750, 1000];
VIEWS = ["Direct", "Backwall Skip", "Sidewall Skip", "Back-Side Skip"];
MODEL = ["NGeo", "Geo"];

times = [];
config = [];
FUNC = [];
xwidth = 60.0e-3;
N = xwidth * 10^3; % + 1;
half_probe_width = 16.0e-3;
xmins = linspace(-xwidth+half_probe_width, -half_probe_width, xwidth*10^3);
xmaxs = linspace(+half_probe_width, xwidth-half_probe_width, xwidth*10^3);

% for i = 48:60%:size(xmins,2)
%     shape = [xmins(i), xmaxs(i), 0.0e-3, 40.0e-3];
%     func = fn_batch_sens_contact_quicker(1.00e-3, 2.0e-3, 500, 3, 0, shape);
%     times = [times; func.times];
%     config = [config; [2, 0, 1.00e-3, 2.0e-3, 500]];
% end

% shape = [xmins(48), xmaxs(48), 0.0e-3, 40.0e-3];
% func = fn_batch_sens_contact_quicker(1.00e-3, 5.0e-3, 500, 2, 0, shape);
shape = [-25e-3, 25e-3, 0.0e-3, 40.0e-3];
fn_batch_sens_contact(1.00e-3, 20e-3, 500, 1, 0, 1);

in_use = monitor_memory_whos;
disp(in_use)

% for i = 1:3
%     for j = 1:2
%         for k = 1:2
%             for view = 1 : 3
%                 for geo = 0 : 1
%                     filename = sprintf('Sens Contact %s %s - p=%.2e pix=%.2e walls=%d.fig', VIEWS(view), MODEL(geo+1), PITCH(i), PIXEL(j), WALLS(k));
%                     if isequal(exist(filename, 'file'), 2)
%                         continue
%                     end
%                     
%                     func = fn_batch_sens_contact(PITCH(i), PIXEL(j), WALLS(k), view, geo);
%                     FUNC = [FUNC; func];
%                     times = [times; func.times];
%                     config = [config; [view, geo, PITCH(i), PIXEL(j), WALLS(k)]];
%                 end
%             end
%         end
%     end
% end

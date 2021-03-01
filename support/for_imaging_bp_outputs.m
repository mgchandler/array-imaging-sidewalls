cd("C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\1 - sensitivity\Figures\Contact\Sensitivity\BP")
dd = dir('*.fig');
filenames = {dd.name};

xwidth = 100.0e-3;
N = 100;
half_probe_width = 1.00e-3 * 32 / 2;
xmins = linspace(-xwidth+half_probe_width, -half_probe_width, N);
xmaxs = linspace(+half_probe_width, xwidth-half_probe_width, N);

openname = filenames{1};
savename = openname(92:end-5);

for f = 1:size(filenames, 2)
    openname = filenames{f};
    savename = sprintf('%s.png', openname(1:end-4));
    s(f).name = savename;
    openfig(openname);
    saveas(gcf, savename);
    close(gcf);
end

count = 0;
less_count = N;
for n = 1:N
    name = sprintf('%.2e,%.2e,%.2e,%.2e', xmins(n), xmaxs(n), 0, 0.04);
    isequal = 0;
    for f = 1:size(s, 2)
        if name == s(f).name
            isequal = 1;
            count = count + 1;
            break
        end
    end
    
    if ~isequal
        fprintf('%d: %s\n', n, name)
        less_count = less_count - 1;
    end
end
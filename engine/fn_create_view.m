function view = fn_create_view(path1, path2, varargin)
% Assembles a view from two paths.
%
% INPUTS:
% - path1 : struct (1, 1)
%       All information of the first path, which will be treated as the
%       forward path in this view. Must be an output from fn_compute_ray
%       function.
% - path2 : struct (1, 1)
%       All information of the first path, which will be treated as the
%       reverse path in this view. Must be an output from fn_compute_ray
%       function.
% - scat_info : OPTIONAL struct (1, 1)
%       Contains information on all scatterers which are being modelled.
%       Scatterer coordinates are in field `image_block`. Must be an output
%       from the fn_scat_info function.
%
% OUTPUTS:
% - view : struct (1, 1)
%       View which collects two paths into one view.

[probe_els, num_scatterers] = size(path1.min_times);

% Unpack paths.
view.path_1 = path1;
view.path_2 = path2;
view.name = sprintf('%s-%s', path1.path_info.name, reverse(path2.path_info.name));

% Initialise variables.
view.min_times = zeros(probe_els ^ 2, num_scatterers);
view.probe_txrx = zeros(probe_els ^ 2, 2);
view.tau = zeros(probe_els ^ 2, path1.scat_size(2), path1.scat_size(1));
view.scatterer_coords = zeros(num_scatterers, 3);

% Assemble view from paths.
for i = 1 : probe_els
    for j = 1 : probe_els
        el = probe_els * (i - 1) + j;
        view.min_times(el, :) = path1.min_times(i, :) + path2.min_times(j, :);
        view.probe_txrx(el, 1) = i;
        view.probe_txrx(el, 2) = j;
        view.tau(el, :, :) = reshape(view.min_times(el, :), [path1.scat_size(2), path1.scat_size(1)]);
    end
end

% If weights have already been calculated, then assemble them all into one.
if isfield(path1, "weights") && isfield(path2, "weights")
    [~, ~, num_freqs] = size(path1.weights.weights);
    view.weights = zeros(probe_els^2, num_scatterers, num_freqs);
    for i = 1 : probe_els
        for j = 1 : probe_els
            el = probe_els * (i - 1) + j;
            view.weights(el, :, :) = ( ...
                path1.weights.weights(i, :, :) .* path2.weights.inv_weights(j, :, :) ...
            );
        end
    end
end

% If scatterer info is provided, then calculate the scattering amplitudes.
if nargin > 2
    scat_info = varargin{1};
    view.scat_amps = fn_scattering_amps(view, path1.freq_array);
end

end
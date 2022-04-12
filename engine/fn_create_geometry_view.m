function view = fn_create_geometry_view(path_info, all_geometries, probe_as_scatterer)
% Computes a single view resulting from signal reflecting from the
% geometry, and as a result does not require two paths in the same way that
% fn_create_view does.
%
% INPUTS:
% - path_info : struct
%       Contains the information from paths, must be output from
%       fn_path_info.
% - all_geometries : struct (num_geoms, 1);
%       All of the geometries present, for checking validity of paths.
% - probe_as_scatterer : struct
%       All of the information about the probe, presented as scatterer.
% - probe_freq : double
%       The frequency of the probe.
%
% OUTPUTS:
% - view
%       A single view which reflects from the geometry.

[probe_els, ~] = size(probe_as_scatterer.x); 

% Compute the ray.
ray = fn_compute_ray(probe_as_scatterer, path_info, all_geometries);

% Initialise this view's arrays.
view.name = path_info.name;
view.min_times = zeros(probe_els ^ 2, 1);
view.probe_txrx = zeros(probe_els ^ 2, 2);
view.valid_path = zeros(probe_els ^ 2, 1);
view.ray = ray;

el = 1;
% Assemble view from ray.
for t_el = 1 : probe_els
    for r_el = 1 : probe_els

        view.min_times(el, 1) = ray.min_times(t_el, r_el);
        view.probe_txrx(el, :) = [t_el, r_el];
        view.valid_path(el) = ray.valid_paths(t_el) .* ray.valid_paths(r_el);

        el = el+1;
    end
end

% % Get the ray weights.
% [~, ~, num_freqs] = size(ray.weights.weights);
% view.weights = zeros(probe_els^2, num_freqs);
% el = 1;
% for ii = 1 : probe_els
%     for jj = 1 : probe_els
%         view.weights(el, :) = ( ...
%             ray.weights.weights(ii, jj, :) * ray.weights.inv_directivity(jj, ii, :) ...
%         );
% 
%     el = el+1;
%     end
% end

end
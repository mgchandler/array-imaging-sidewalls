function view = fn_create_view(path1, path2)
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
%
% OUTPUTS:
% - view : struct (1, 1)
%       View which collects two paths into one view.

[probe_els, num_scatterers] = size(path1.min_times);
if ~strcmp(path1.scat_info.type, "image")
    fmc_mask = path1.scat_info.fmc_mask;
end
% Memory limit above which to switch to single precision to halve memory
% requirement of min_times, ray_weights and scat_amps in view.
mb_limit = 100;

% Unpack paths.
view.path_1 = path1;
view.path_2 = path2;
view.name = sprintf('%s - %s', path1.path_info.name, path2.path_info.rev_name);

% Initialise variables.
if 8*probe_els^2*num_scatterers / 1024^2 > mb_limit
    view.min_times = single(zeros(probe_els ^ 2, num_scatterers));
    view.mask = single(ones(probe_els ^ 2, num_scatterers));
else
    view.min_times = zeros(probe_els ^ 2, num_scatterers);
    view.mask = ones(probe_els ^ 2, num_scatterers);
end
view.probe_txrx = int8(zeros(probe_els ^ 2, 2));
view.valid_path = logical(zeros(probe_els ^ 2, num_scatterers));
% view.mask = ones(probe_els ^ 2, num_scatterers);
% view.tx_angles = zeros(probe_els, num_scatterers);
% view.rx_angles = zeros(probe_els, num_scatterers);
% view.scat_inc_angles = zeros(probe_els, num_scatterers);
% view.scat_out_angles = zeros(probe_els, num_scatterers);

% Assemble view from paths.
for ii = 1 : probe_els
    for jj = 1 : probe_els
        el = probe_els * (ii - 1) + jj;
        view.min_times(el, :) = path1.min_times(ii, :) + path2.min_times(jj, :);
        view.probe_txrx(el, 1) = ii;
        view.probe_txrx(el, 2) = jj;
        % If path is not valid, it will be 0. If 1, path is valid. View
        % path will be invalid if either path is invalid, so multiply.
        view.valid_path(el, :) = and(path1.valid_paths(ii, :), path2.valid_paths(jj, :));
        
%         if isfield(path2, "weights")
%             if ii == probe_els
%                 view.rx_angles(jj, :) = path2.weights.out_theta(jj, :, 1, 2);
%                 view.scat_out_angles(jj, :) = path2.weights.inv_out_theta(jj, :, 1, 2);
%             end
%         end
    end
%     if isfield(path1, "weights")
%         view.tx_angles(ii, :) = path1.weights.out_theta(ii, :, 1, 2);
%         view.scat_inc_angles(ii, :) = path1.weights.inc_theta(ii, :, end, 2);
%     end
end

% If weights have already been calculated, then assemble them all into one.
if isfield(path1, "weights") && isfield(path2, "weights")
    [~, ~, num_freqs] = size(path1.weights.weights);
    view.weights = zeros(probe_els^2, num_scatterers, num_freqs);
    if 8*probe_els^2*num_scatterers / 1024^2 > mb_limit
        view.weights = single(view.weights);
    end
    el = 0;
    for ii = 1 : probe_els
        for jj = 1 : probe_els
            el = el+1;
            if isstruct(fmc_mask)
                view.mask(el, :) = interp1(fmc_mask.time, fmc_mask.data(:, el), view.min_times(el, :));
            end
            view.weights(el, :, :) = ( ...
                path1.weights.weights(ii, :, :) .* path2.weights.inv_weights(jj, :, :) .* view.mask(el, :) ...
            );
        end
    end
end

% If scatterer info is provided, then calculate the scattering amplitudes.
if isfield(path1, 'freq_array')
    view.scat_amps = conj(fn_scattering_amps(view, path1.freq_array));
    if 8*probe_els^2*num_scatterers / 1024^2 > mb_limit
        view.scat_amps = single(view.scat_amps);
    end
end

end
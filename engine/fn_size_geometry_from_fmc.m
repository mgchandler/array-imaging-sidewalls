function geometry = fn_size_geometry_from_fmc(exp_data, varargin)
% Works out the exact coordinates of the L-shaped geometry from supplied
% FMC data. Assumes that exp_data is saved by Brain (i.e. contains
% attributes "time", "time_data", "array", "material", "tx", "rx")
% exp_data expected output from Brain read from array.
%
% Inputs:
%   - exp_data : contents of file saved by Brain in FMC.
%   - z_depth  : Approximate depth of the back wall below the array (mm)
%   - x_depth  : Approximate distance of the side wall from the centre of
%                the array (mm)

boxsize = 5e-3;
n_pix = 50;

%% Set up settings to be reused throughout the test.
if nargin == 1
    z_depth = abs(input("Approximately how deep is the backwall from the array? (mm)\n>> ")) * 1e-3;
    x_depth = abs(input("Approximately how far is the sidewall from the array? (mm)\n>> ")) * 1e-3;
else
    z_depth = varargin{1};
    x_depth = varargin{2};
    if nargin > 3
        is_r = varargin{3};
        if nargin > 4
            savename = varargin{4};
        end
    end
end
bw_settings.probe.freq = exp_data.array.centre_freq;
bw_settings.probe.pitch = exp_data.array.el_xc(2) - exp_data.array.el_xc(1);
bw_settings.probe.num_els = max(exp_data.tx);
bw_settings.model.pixel = boxsize/n_pix;
bw_settings.model.plot_it = false;
freq = [0:length(exp_data.time)-1] / exp_data.time(end);

%% Start by imaging near to the back wall
bw_settings.mesh.geom.x = {min(exp_data.array.el_xc) - boxsize, min(exp_data.array.el_xc) - boxsize, max(exp_data.array.el_xc) + boxsize, max(exp_data.array.el_xc) + boxsize};
bw_settings.mesh.geom.z = {0, z_depth*1.5, z_depth*1.5, 0};
bw_settings.model.image_range = [-boxsize/2, boxsize/2, z_depth - boxsize/2, z_depth + boxsize/2];
bw_settings.model.max_no_reflections = 0;
bw_settings.data.time = exp_data.time - 5e-7;
HMC_time_data = ifft(2 * fn_hanning(length(exp_data.time), bw_settings.probe.freq/max(freq), bw_settings.probe.freq/max(freq)) .* fft(exp_data.time_data));
if size(HMC_time_data, 2) ~= max(exp_data.tx)^2
    [FMC_time_data, ~, ~] = fn_swap_FMC_HMC(exp_data.time_data, exp_data.tx, exp_data.rx);
else
    FMC_time_data = HMC_time_data;
end
bw_settings.data.data = FMC_time_data;

model_options = fn_default_model_options(bw_settings);
[bw_Ims, ~, ~] = fn_tfm(model_options);

signals_db  = zeros(size(bw_Ims(1).image, 1), 1);
for im = 1:size(bw_Ims, 1)
    signals_db  = signals_db  + sum(bw_Ims(im).db_image, 2);
end
bw_depth = bw_Ims(1).z(signals_db == max(signals_db));

%% Next get the side wall.
sw_settings = bw_settings;
sw_settings.mesh = rmfield(sw_settings.mesh, 'geom');
sw_settings.model.max_no_reflections = 1;
sw_settings.model.wall_for_imaging = 'B1';
sw_settings.mesh.geom.z = {0, bw_depth, bw_depth, 0};
sw_settings.model.pixel = 2*boxsize/n_pix;
if exist('is_r', 'var') % Only do the correct one
    if is_r
        sw_settings.model.image_range = [x_depth - boxsize, x_depth + boxsize, boxsize, bw_depth - boxsize];
        sw_settings.mesh.geom.x = {min(exp_data.array.el_xc) - boxsize, min(exp_data.array.el_xc) - boxsize, x_depth*1.5, x_depth*1.5};
    else
        sw_settings.model.image_range = [-x_depth - boxsize, -x_depth + boxsize, boxsize, bw_depth - boxsize];
        sw_settings.mesh.geom.x = {-x_depth*1.5, -x_depth*1.5, max(exp_data.array.el_xc) + boxsize, max(exp_data.array.el_xc) + boxsize};
    end
    model_options = fn_default_model_options(sw_settings);
    [sw_Ims, ~, ~] = fn_tfm(model_options);
else % Do both and work out which side to use
    sw_settings.model.image_range = [x_depth - boxsize, x_depth + boxsize, boxsize, bw_depth - boxsize];
    sw_settings.mesh.geom.x = {min(exp_data.array.el_xc) - boxsize, min(exp_data.array.el_xc) - boxsize, x_depth*1.5, x_depth*1.5};
    model_options = fn_default_model_options(sw_settings);
    [sw_ImsR, ~, ~] = fn_tfm(model_options);
    
    sw_settings.model.image_range = [-x_depth - boxsize, -x_depth + boxsize, boxsize, bw_depth - boxsize];
    sw_settings.mesh.geom.x = {-x_depth*1.5, -x_depth*1.5, max(exp_data.array.el_xc) + boxsize, max(exp_data.array.el_xc) + boxsize};
    model_options = fn_default_model_options(sw_settings);
    [sw_ImsL, ~, ~] = fn_tfm(model_options);

    max_ = 0;
    is_r = true;
    for im = 1:size(sw_ImsR, 1)
        [max_im, l_or_r] = max([max(abs(sw_ImsL(im).image(:))), max(abs(sw_ImsR(im).image(:)))]);
        if max_im > max_
            max_ = max_im;
            is_r = logical(l_or_r - 1);
        end
    end

    if is_r
        sw_Ims = sw_ImsR;
    else
        sw_Ims = sw_ImsL;
    end
end

signals_db  = zeros(1, size(sw_Ims(1).image, 2));
for im = 1:size(sw_Ims, 1)
    % Only do up to half-skip.
    if count(sw_Ims(im).name, "L") + count(sw_Ims(im).name, "T") < 4
        signals_db  = signals_db + sum(sw_Ims(im).db_image, 1);
    end
end
sw_depth = sw_Ims(1).x(signals_db == max(signals_db));

%% Finally get nlos corner
nlos_settings = sw_settings;
nlos_settings.mesh = rmfield(nlos_settings.mesh, 'geom');
if is_r
    nlos_settings.model.wall_for_imaging = "S1";
    nlos_settings.mesh.geom.x = {sw_depth, sw_depth, 0, 0, min(exp_data.array.el_xc) - boxsize, min(exp_data.array.el_xc) - boxsize};
    nlos_settings.mesh.geom.z = {0, bw_depth + 25e-3, bw_depth + 25e-3, bw_depth, bw_depth, 0};
    nlos_settings.model.image_range = [sw_depth - 20e-3 - boxsize, sw_depth - 20e-3 + boxsize, bw_depth + 20e-3 - boxsize, bw_depth + 20e-3 + boxsize];
else
    nlos_settings.model.wall_for_imaging = "S3";
    nlos_settings.mesh.geom.x = {max(exp_data.array.el_xc) + boxsize, max(exp_data.array.el_xc) + boxsize, 0, 0, sw_depth, sw_depth};
    nlos_settings.mesh.geom.z = {0, bw_depth, bw_depth, bw_depth + 25e-3, bw_depth + 25e-3, 0};
    nlos_settings.model.image_range = [sw_depth + 20e-3 - boxsize, sw_depth + 20e-3 + boxsize, bw_depth + 20e-3 - boxsize, bw_depth + 20e-3 + boxsize];
end
model_options = fn_default_model_options(nlos_settings);
[nlos_Ims, ~, ~] = fn_tfm(model_options);

fuse = 0;
for ii = 1:size(nlos_Ims, 1)
    if count(nlos_Ims(im).name, "L") + count(nlos_Ims(im).name, "T") == 4
        fuse = fuse + nlos_Ims(ii).db_image;
    end
end
[~, ind] = max(fuse(:));
[zi, xi] = ind2sub(size(nlos_Ims(1).db_image), ind);
nlos_x = nlos_Ims(1).x(xi); nlos_z = nlos_Ims(1).z(zi);

if is_r
    geom_corners = [sw_depth, sw_depth, nlos_x, nlos_x, min(exp_data.array.el_xc) - boxsize, min(exp_data.array.el_xc) - boxsize;
                    zeros(1,6);
                    0, nlos_z, nlos_z, bw_depth, bw_depth, 0].';
else
    geom_corners = [min(exp_data.array.el_xc) - boxsize, min(exp_data.array.el_xc) - boxsize, nlos_x, nlos_x, sw_depth, sw_depth;
                    zeros(1,6);
                    0, bw_depth, bw_depth, nlos_z, nlos_z, 0].';
end

geometry = fn_make_geometry(true, model_options.mesh.geom.profile_list, model_options.mesh.geom.n_pts, ...
    geom_corners ...
);

save(savename, 'bw_Ims', 'sw_Ims', 'nlos_Ims', 'geometry')

end
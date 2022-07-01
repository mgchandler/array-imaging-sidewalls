clear
clc
% close all

% profile clear
% profile -memory on -historysize 25000000

% Are we using data generated from BP, or are we running
% array-imaging-sidewalls?
is_bp_data = true;
% If we're running a-i-s, are we using tfm or sens to get our signal
% values? N.B. If is_bp_data = 1, then this logical is not used.
is_tfm = true;
% Are we modelling with geometry or not?
is_geom = false;
% Use analytical wave velocities?
is_book_velocity = false;

image_block = [[0.0e-3, 0.0, 17.5e-3]; ...
     [16.25e-3, 0.0, 17.5e-3]; ...
     [32.5e-3, 0.0, 5.625e-3]; ...
     [32.5e-3, 0.0, 16.875e-3]; ...
     [32.5e-3, 0.0, 28.125e-3]; ...
     [32.5e-3, 0.0, 39.375e-3]];

% Measured speeds from FE model.
v_L = 6305.350479;%6317.0122248907810;
v_S = 3152.675239;%3110.2818131859126;

if is_bp_data
    cd("C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\Abaqus\AbaqusInputFileGeneration - Output\v9\Output\FMC Data\FE - L\old_20_for_comparison")
else
    cd("C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\Abaqus\AbaqusInputFileGeneration - Output\v9\Output\FMC Data\FE - L\old_20_for_comparison")
end

npw = [15:5:60];
npw = 20;

Views_im = 0;
maxes = zeros(size(npw, 2), 6, 21);
db_maxes = zeros(size(npw, 2), 6, 21);

if is_geom
    geom_str = 'geom';
else
    geom_str = 'nogeom';
end

scats_to_run = 4:6;

if or(is_bp_data, is_tfm)
    for ii = scats_to_run
        yaml_options = yaml.loadFile(sprintf("L_%d.yml", ii));
        yaml_options.material.couplant_v = 340.0;
        yaml_options.material.couplant_density = 1.2;
        for kk = 1:size(yaml_options.mesh.geom.z, 2)
            yaml_options.mesh.geom.z{kk} = -yaml_options.mesh.geom.z{kk};
        end
        yaml_options.model.boxsize = 1.75e-3;
        yaml_options.model.interp_method = 'lanczos';
        yaml_options.model.pixel = .1e-3;
        yaml_options.model.model_geom = is_geom;
        yaml_options.model.wall_for_imaging = "S1";
        yaml_options.probe.angle = 0.0;
        yaml_options.probe.standoff = 0.0;
%         yaml_options.mesh.n_per_wl = npw;
        sdh_rad = yaml_options.mesh.scat.r{1};

        if is_bp_data
            yaml_options.model.savepath = "C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\Abaqus\AbaqusInputFileGeneration - Output\v9\Output\FMC Data\FE - L\old_20_for_comparison";
        else
            yaml_options.model.savepath = "C:\Users\mc16535\OneDrive - University of Bristol\Documents\Postgrad\Coding\Abaqus\AbaqusInputFileGeneration - Output\v9\Output\FMC Data\FE - L\old_20_for_comparison";
        end
        
        filename = sprintf('L_%dnpw_%d', npw, ii);
        disp(filename)
        if is_bp_data
            if ~is_geom
                load(sprintf('L_%d_bl_BP.mat', ii));
                Bl_data = data;
            end
            load(sprintf('L_%d_BP.mat', ii))
            if ~is_book_velocity
                v_L = fn_speed_from_fmc(time(:, 1), data, reshape(repmat(1:32, 32, 1), 1, 1024), repmat(1:32, 1, 32), 25.e-3);
                v_S = .5 * v_L;
            else
                v_L = 6317.0122248907810;
                v_S = 3110.2818131859126;
            end
            if ~is_geom
                data = data - Bl_data(1:size(data, 1), 1:size(data, 2));
            end
        else
            v_L = 6317.0122248907810;
            v_S = 3110.2818131859126;
        end
        scat_coords = image_block(ii, :);

        yaml_options.mesh.scat = fn_scat_info( ...
            "sdh", ...
            scat_coords(1), ...
            scat_coords(2), ...
            scat_coords(3), ...
            sdh_rad, ...
            v_L/5e6, ...
            v_S/5e6, ...
            deg2rad(0), ...
            'ang_pts_over_2pi', 120 ...
        );

        image_box = 4e-3;
        yaml_options.model.image_range = [scat_coords(1)-image_box, scat_coords(1)+image_box, scat_coords(3)-image_box, scat_coords(3)+image_box];

        xsize = yaml_options.model.image_range(2) - yaml_options.model.image_range(1);
        zsize = yaml_options.model.image_range(4) - yaml_options.model.image_range(3);
        xpts = round(xsize / yaml_options.model.pixel);
        zpts = round(zsize / yaml_options.model.pixel);
        x = linspace(yaml_options.model.image_range(1), yaml_options.model.image_range(2), xpts+1);
        z = linspace(yaml_options.model.image_range(3), yaml_options.model.image_range(4), zpts+1);
        [X, Z] = meshgrid(x, z);
            
        if is_bp_data
            yaml_options.data.time = squeeze(time(:, 1));
            yaml_options.data.data = data;
            yaml_options.model.savename = strcat(filename, '_', geom_str, '_measV');

            filename = yaml_options.model.savename;
            yaml_options.model.max_no_reflections = 1;
            model_options = fn_default_model_options(yaml_options);
            
            if and(isstruct(Views_im), is_book_velocity)
                [Ims, Views_im] = fn_tfm(model_options, Views_im);
            else
                [Ims, Views_im] = fn_tfm(model_options);
            end

%             xsize = model_options.model.image_range(2) - model_options.model.image_range(1);
%             zsize = model_options.model.image_range(4) - model_options.model.image_range(3);
%             xpts = round(xsize / model_options.model.pixel);
%             zpts = round(zsize / model_options.model.pixel);
%             x = linspace(model_options.model.image_range(1), model_options.model.image_range(2), xpts+1);
%             z = linspace(model_options.model.image_range(3), model_options.model.image_range(4), zpts+1);
%             [X, Z] = meshgrid(x, z);
            
%             load(strcat(filename, '.mat'))

            box_mask = and(and(X>scat_coords(1)-model_options.model.boxsize/2, X<scat_coords(1)+model_options.model.boxsize/2), ...
                           and(Z>scat_coords(3)-model_options.model.boxsize/2, Z<scat_coords(3)+model_options.model.boxsize/2));

            for im = 1:21
                maxes(1, ii, im) = max(Ims(im).image(box_mask));
            end
            db_maxes(1, ii, :) = 20 * log10(abs(maxes(1, ii, :)) ./ abs(maxes(1, ii, 12))); 
            disp(" ")
        else
            try
                yaml_options = rmfield(yaml_options, 'FMC_data');
            end
            
            which = 1;

            yaml_options.model.savename = strcat(filename, '_', geom_str, '_measV');
            yaml_options.model.npw = npw;

            filename = yaml_options.model.savename;
            model_options = fn_default_model_options(yaml_options);
            
            if and(isstruct(Views_im), is_book_velocity)
                [Ims, Views_im] = fn_tfm(model_options, Views_im);
            else
                [Ims, Views_im] = fn_tfm(model_options);
            end

%             xsize = model_options.model.image_range(2) - model_options.model.image_range(1);
%             zsize = model_options.model.image_range(4) - model_options.model.image_range(3);
%             xpts = round(xsize / model_options.model.pixel);
%             zpts = round(zsize / model_options.model.pixel);
%             x = linspace(model_options.model.image_range(1), model_options.model.image_range(2), xpts+1);
%             z = linspace(model_options.model.image_range(3), model_options.model.image_range(4), zpts+1);
%             [X, Z] = meshgrid(x, z);
            
            for which_npw = 1:size(npw, 2)
%                 load(strcat(filename, '.mat'))

                box_mask = and(and(X>scat_coords(1)-model_options.model.boxsize/2, X<scat_coords(1)+model_options.model.boxsize/2), ...
                               and(Z>scat_coords(3)-model_options.model.boxsize/2, Z<scat_coords(3)+model_options.model.boxsize/2));

                for im = 1:21
                    maxes(which_npw, ii, im) = max(Ims(im).image(box_mask));
                end
                db_maxes(which_npw, ii, :) = 20 * log10(abs(maxes(which_npw, ii, :)) ./ abs(maxes(which_npw, ii, 12))); 
                disp(" ")
            end
        end
    end

else

    filename = strcat('L_sens_', sprintf('%dnpw', npw));
    
    yaml_options.model.image_locs = image_block;

    yaml_options.mesh.scat = fn_scat_info( ...
        "sdh", ...
        0, ...
        0, ...
        0, ...
        sdh_rad, ...
        v_L/5e6, ...
        v_S/5e6, ...
        deg2rad(0), ...
        'ang_pts_over_2pi', 120 ...
    );
    yaml_options.model.savename = filename;
    model_options = fn_default_model_options(yaml_options);
    fn_sens(model_options)

    load(strcat(filename, '.mat'))
    
    for which_npw = 1:size(npw,2)
        if size(model_options.model.image_locs, 2) == 3
            for im = 1:21
                maxes(which_npw, :, im) = Sens(im).image;
            end
        else
            xsize = model_options.model.image_range(2) - model_options.model.image_range(1);
            zsize = model_options.model.image_range(4) - model_options.model.image_range(3);
            xpts = round(xsize / model_options.model.pixel);
            zpts = round(zsize / model_options.model.pixel);
            x = linspace(model_options.model.image_range(1), model_options.model.image_range(2), xpts+1);
            z = linspace(model_options.model.image_range(3), model_options.model.image_range(4), zpts+1);
            [X, Z] = meshgrid(x, z);
            for ii = 1:size(image_block, 1)
                for im = 1:21
                    maxes(which_npw, ii, im) = interp2(X, Z, Sens(im).image, image_block(ii, 1), image_block(ii, 3), 'cubic');
                end
            end
        end

        for scat = scats_to_run
            db_maxes(which_npw, scat, :) = 20 * log10(abs(maxes(which_npw, scat, :)) ./ abs(maxes(which_npw, scat, 12))); 
        end
    end

end

real_maxes = real(maxes);
imag_maxes = imag(maxes);
abs_maxes  = squeeze(abs(maxes));
db_maxes = squeeze(db_maxes);

% profile report
VIEWS = 2;
Number_of_ims = 21;
plot_x = 3;
plot_z = 7;
im_width = 450;
im_height = 1050;
UC = 10^3;
db_range_for_output = 40;

PITCH = model_config.PITCH;
PIXEL = model_config.PIXEL;
WALLS = model_config.WALLS;
VIEWS = model_config.VIEWS;
GEOM = model_config.GEOM;
SETUP = model_config.SETUP;

probe_els = model_options.probe_els;
probe_angle = model_options.probe_angle;
probe_standoff = model_options.probe_standoff;
probe_frequency = model_options.probe_frequency;
el_length = model_options.el_length + PITCH;
xmin = model_options.geom_shape.xmin;
xmax = model_options.geom_shape.xmax;
zmin = model_options.geom_shape.zmin;
zmax = model_options.geom_shape.zmax;
couplant_speed = model_options.material_params.couplant_speed;
couplant_density = model_options.material_params.couplant_density;
solid_long_speed = model_options.material_params.solid_long_speed;
solid_shear_speed = model_options.material_params.solid_shear_speed;
solid_density = model_options.material_params.solid_density;
boxsize = PIXEL * model_options.boxsize;
image_block_info = model_options.scat_info;
savepath = model_options.savepath;
savename = model_options.savename;

front_wall_pos = zmin - 1.0e-5;
front_wall_pixels = WALLS;
back_wall_pos = zmax + 1.0e-5;
back_wall_pixels = WALLS;
side_wall_pos = xmax + 1.0e-5;
side_wall_pixels = WALLS;

rot_matrix = [cos(probe_angle) 0 sin(probe_angle); 0 1 0; -sin(probe_angle) 0 cos(probe_angle)];

probe_coords = zeros(3, probe_els);
probe_coords(1, :) = linspace(0, (probe_els - 1) * el_length, probe_els);
probe_coords(1, :) = probe_coords(1, :) - mean(probe_coords(1, :));
probe_coords = probe_coords.' * rot_matrix;
probe_coords(:, 3) = probe_coords(:, 3) - probe_standoff;

front_wall = zeros(front_wall_pixels, 3);
front_wall(:, 1) = linspace(xmin, xmax, front_wall_pixels);
front_wall(:, 3) = linspace(front_wall_pos, front_wall_pos, front_wall_pixels);

back_wall = zeros(back_wall_pixels, 3);
back_wall(:, 1) = linspace(xmin, xmax, front_wall_pixels);
back_wall(:, 3) = linspace(back_wall_pos, back_wall_pos, back_wall_pixels);

side_wall = zeros(side_wall_pixels, 3);
side_wall(:, 1) = linspace(side_wall_pos, side_wall_pos, side_wall_pixels);
side_wall(:, 3) = linspace(zmin, zmax, side_wall_pixels);



xsize = xmax - xmin;
zsize = zmax - zmin;

xpts = round(xsize / PIXEL);
zpts = round(zsize / PIXEL);

% Sort into same order as arim for easy comparison.
if VIEWS == 1
    View_names = ["L-L", "L-T", "T-T"];
elseif VIEWS == 2
    View_names = ["L-L", "L-T", "T-T", "L-LBL", "T-LBL", "L-TBL", "T-TBL", ...
                  "L-LBT", "T-LBT", "L-TBT", "T-TBT", "LBL-LBL", "LBL-LBT", "LBL-TBL", ...
                  "LBL-TBT", "LBT-LBT", "LBT-TBL", "LBT-TBT", "TBL-LBT", "TBL-TBT", "TBT-TBT"];
elseif VIEWS == 3
    View_names = ["L-L", "L-T", "T-T", "L-LSL", "T-LSL", "L-TSL", "T-TSL", ...
                  "L-LST", "T-LST", "L-TST", "T-TST", "LSL-LSL", "LSL-LST", "LSL-TSL", ...
                  "LSL-TST", "LST-LST", "LST-TSL", "LST-TST", "TSL-LST", "TSL-TST", "TST-TST"];
elseif VIEWS == 4
    View_names = ["L-L", "L-T", "T-T", "LBL-T", "LBL-T", "LBT-L", "LBT-T", ...
                  "TBL-L", "TBL-T", "TBT-L", "TBT-T", "LBL-LBL", "LBL-LBT", "LBL-TBL", ...
                  "LBL-TBT", "LBT-LBT", "LBT-TBL", "LBT-TBT", "TBL-LBT", "TBL-TBT", "TBT-TBT", ...
                  "L-LSL", "T-LSL", "L-TSL", "T-TSL", ...
                  "L-LST", "T-LST", "L-TST", "T-TST", "LSL-LSL", "LSL-LST", "LSL-TSL", ...
                  "LSL-TST", "LST-LST", "LST-TSL", "LST-TST", "TSL-LST", "TSL-TST", "TST-TST"];
end

sort_idx = zeros(Number_of_ims, 1);
view_names = struct2cell(Views);
view_names = view_names(3, :);
for view = 1:Number_of_ims
    sort_idx(view) = find(strcmp([view_names{:}], View_names(view)));
end

% Plot.
im_x = linspace(xmin, xmax, xpts+1);
im_z = linspace(zmin, zmax, zpts+1);
fig = figure(1);
ax = repmat(subplot(plot_z, plot_x, 1), Number_of_ims, 1);
if image_block_info.type == "crack"
    sgtitle(sprintf('Sens %.2f Crack - %.2f deg', image_block_info.crack_length, rad2deg(image_block_info.angle)))
end
for im = 1:Number_of_ims
    ax(im) = subplot(plot_z, plot_x, im);
    imagesc(im_x*UC, im_z*UC, Sens(sort_idx(im)).db_image);
    hold on
    title(Sens(sort_idx(im)).name)
    caxis([-db_range_for_output, 0])
    plot(probe_coords(:, 1)*UC, probe_coords(:, 3)*UC, 'go');
    plot(front_wall(:, 1)*UC, front_wall(:, 3)*UC, 'r');
    plot(back_wall(:, 1)*UC, back_wall(:, 3)*UC, 'r');
    plot(side_wall(:, 1)*UC, side_wall(:, 3)*UC, 'r');
    axis equal; axis tight;
end

% Resize image and move for colorbar
set(fig, 'Position', [20, 20, im_width, im_height])
posa = cell2mat(get(ax, 'Position'));
h = colorbar;
set(ax(im), 'Position', posa(im, :))
set(ax, 'units', 'pix')
set(h, 'units', 'pix')
posf = get(fig, 'Position'); % gives x left, y bottom, width, height
set(fig, 'Position',  [posf(1:2) posf(3)*1.1 posf(4)])
hpos = h.Position;
posa = cell2mat(get(ax, 'Position'));
pos1 = posa(1, :);
set(h, 'Position', [hpos(1)+10, hpos(2), hpos(3)*2, pos1(2)-hpos(2)+pos1(4)])

h.Label.String = 'dB';
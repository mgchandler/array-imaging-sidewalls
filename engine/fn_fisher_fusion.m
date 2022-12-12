function Fused = fn_fisher_fusion(model_options, Views)
% Performs fusion using Fisher's method for Views which have previously
% been obtained, containing probabilities that the intensity at a pixel
% location comes from an artefact (H0: intensity is from an artefact). Note
% that grain noise is assumed non-existent thus σ=1.
%
% INPUTS:
% - model_options : struct
%       Output of fn_default_model_options(); presently only used for save
%       location and save name.
% - Views : Nx1 struct
%       Modified version of Views from fn_make_views(). Assumes that
%       db_image field contains probabilities, and image field logicals
%       expressing whether to include the pixel in the fusion process
%       (should either be 1 or ndarray of logicals of equal size to
%       db_image).

%% Fuse
number_of_ims = size(Views, 1);
Fused = Views(1);
Fused.name = "Fused TFM";
Fused.q = 0;
dof = zeros(size(Views(1).image));
for im = 1:number_of_ims
    if or(any(Views(im).db_image < 0), any(Views(im).db_image > 1))
        error("fn_fisher_fusion: Views do not contain probabilities.")
    end
    Views(im).p = Views(im).db_image;
    % Calculate q: input parameter to chi2 cdf.
    q_im = 2*log(Views(im).p);
    q_im(~Views(im).image) = 0;
    Fused.q = Fused.q - q_im;
    dof = dof + Views(im).image;
end
% Fused.p = 1 - chi2cdf(Fused.q, 2*size(Views, 1));
Fused.p = 1 - chi2cdf(Fused.q, 2*dof);

Fused.db_image = Fused.p;
Fused.plotExtras = Views(1).plotExtras;

% for im = 1:number_of_ims
%     Fused.image = Fused.image - 2*log(Views(im).db_image);
% end
% Fused.db_image = chi2cdf(Fused.image, 2*size(Views, 1));

p_threshold = .005;
fn_image_from_mat(Fused)
grp = get(get(gcf, 'Children'), 'Children');
grp(2).CLim = [0, p_threshold];
grp(1).Label.String = "p";

% p_threshold = .95;
% fn_image_from_mat(Fused)
% grp = get(get(gcf, 'Children'), 'Children');
% grp(2).CLim = [p_threshold, 1];
% grp(1).Label.String = "p";

if model_options.model.savepath ~= ""
%     cd(model_options.model.savepath)
    savename = sprintf("%s masked", model_options.model.savename);
    filename_fig = fullfile(model_options.model.savepath, sprintf('%s.fig', savename));
    filename_mat = fullfile(model_options.model.savepath, sprintf('%s.mat', savename));
    savefig(filename_fig)
    filename_fig = fullfile(model_options.model.savepath, sprintf('%s.png', savename));
    saveas(gcf, filename_fig)
    save(filename_mat, 'Fused')
end

end
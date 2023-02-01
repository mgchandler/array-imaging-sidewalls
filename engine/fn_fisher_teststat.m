function Fused = fn_fisher_teststat(model_options, Ims, Sigma)
% Performs Fisher's fusion method for otherwise unprocessed TFMs. Assumes
% that pixel-wise amplitudes are well described by Rician distributions
% which have already been obtained for this TFM's focal law. Note that this
% means that grain noise is assumed to be non-existant.
%
% db_image field is used throughout all structs to represent their values -
% they are not neccessarily images in dB scale. This is done to enable use
% of fn_plot_image_from_mat() without having to swap fields around. As this
% is a preferential thing, more descriptive fields could be used with the
% change to db_image field right before plotting.
%
% INPUTS:
% - model_options : struct
%       Options in form output from fn_default_model_options(). Currently
%       only used for saving location + name.
% - Ims : struct (N, 1)
%       Contains TFMs to be processed.
% - Nu : struct (N, 1)
%       View-wise non-centrality parameter of the Rician distributions.
%       Expected to be of the same form as Ims, however db_image field
%       contains the values for ν.
% - Sigma : struct (N, 1)
%       View-wise scale parameter of the Rician distributions. Expected to
%       be of the same form as Ims, however db_image field contains the
%       values for σ.
% - alpha : double
%       p-level required for rejection of the null hypothesis. This will
%       only be used for plotting. (Default α = 0.005)
% - GoF : struct (N, 1)
%       View-wise evaluation of how well Rician distributions describe
%       intensities. Should be of the same form as Ims, however db_image
%       field contains the probabilities that the null-hypothesis is true.
%       If supplied, the same α used for fusion result will be used here.
%       (Default none supplied).
%
% OUTPUTS:
% - Fused : struct (1, 1)
%       Resulting fused image containing the probability that each pixel is
%       well described by Rician distributions in every view. This means    
%       that for p(r) < α indicates that the pixel at location r is not 
%       well described and might be defective.

number_of_ims = size(Ims, 1);
if size(Sigma, 1) ~= number_of_ims
    error("fn_fisher_fusion: Nu and Sigma must have the same number of views as TFMs supplied")
end

db_range = 40;

Fused = Ims(1);
Fused.name = "Fisher Fusion Test Statistic";
Fused.db_image = 0;
Fused.image = 0;

%% Compute test parameter
for im = 1:number_of_ims
    Fused.image = Fused.image + (Ims(im).image ./ Sigma(im).db_image).^2;
end
Fused.image = sqrt(Fused.image);
Fused.db_image = abs(Fused.image);
%% Save Fused (phased array scale)          
fn_image_from_mat(Fused)
grp = get(get(gcf, 'Children'), 'Children');
grp(2).CLim = [0, 60];
grp(2).XLim = [Ims(1).x(1)*1e3, Ims(1).x(end)*1e3];
grp(2).YLim = [Ims(1).z(1)*1e3, Ims(1).z(end)*1e3];
grp(1).Label.String = "(array units)";
if ~strcmp(model_options.model.savepath, "")
    savefig(fullfile(model_options.model.savepath, sprintf("%s Fisher fusion Test Statistic", model_options.model.savename)))
    saveas(gcf, fullfile(model_options.model.savepath, sprintf("%s Fisher fusion Test Statistic.png", model_options.model.savename)))
    save(fullfile(model_options.model.savepath, sprintf("%s Fisher fusion Test Statistic.mat", model_options.model.savename)), "Fused")
end

max_ = max(abs(Fused.image(:)));
Fused.db_image = 20 * log10(abs(Fused.image) ./ max_); 

%% Save Fused image             
fn_image_from_mat(Fused)
grp = get(get(gcf, 'Children'), 'Children');
grp(2).CLim = [-db_range, 0];
grp(2).XLim = [Ims(1).x(1)*1e3, Ims(1).x(end)*1e3];
grp(2).YLim = [Ims(1).z(1)*1e3, Ims(1).z(end)*1e3];
grp(1).Label.String = "dB";
if ~strcmp(model_options.model.savepath, "")
    savefig(fullfile(model_options.model.savepath, sprintf("%s Fisher fusion Test Statistic (dB)", model_options.model.savename)))
    saveas(gcf, fullfile(model_options.model.savepath, sprintf("%s Fisher fusion Test Statistic (dB).png", model_options.model.savename)))
    save(fullfile(model_options.model.savepath, sprintf("%s Fisher fusion Test Statistic (dB).mat", model_options.model.savename)), "Fused")
end

end
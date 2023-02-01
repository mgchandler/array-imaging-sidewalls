function Fused = fn_fisher_fusion(model_options, Ims, Nu, Sigma, varargin)
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
if or(size(Nu, 1) ~= number_of_ims, size(Sigma, 1) ~= number_of_ims)
    error("fn_fisher_fusion: Nu and Sigma must have the same number of views as TFMs supplied")
end

alpha = 0.005;
GoF = Ims;
for im = 1:number_of_ims
    GoF(im).db_image = ones(size(Ims(im).db_image));
end
if nargin > 4
    for idx = 1:nargin-4
        arg = varargin{idx};
        if isstruct(arg)
            GoF = arg;
        elseif isdouble(arg)
            alpha = arg;
        end
    end
end

Probs = Ims;
Probs = rmfield(Probs, ["db_image", "image"]);

Fused = Ims(1);
Fused.name = "Fisher Fusion";
Fused.db_image = 0;
Fused.image = 0;
dof = zeros(size(Ims(1).db_image));

%% Compute probabilities
for im = 1:number_of_ims
%     GoF(im).db_image(Ims(im).db_image < -40) = 0; % Neglect small amplitudes
    % From RicianDistribution.cdffunc: ncx2cdf((x./sigma).^2, 2, (s./sigma).^2)
    Probs(im).db_image = 1 - ncx2cdf((abs(Ims(im).image)./Sigma(im).db_image).^2, 2, (Nu(im).db_image./Sigma(im).db_image).^2);
    Probs(im).chi = GoF(im).db_image > 0.05;
    %% Fuse
    q = 2*log(Probs(im).db_image);
    q(~Probs(im).chi) = 0;
    Fused.image = Fused.image - q;
    dof = dof + Probs(im).chi;
end
Fused.db_image = 1 - chi2cdf(Fused.image, 2*dof);

%% Save View-wise probabilities and Fused image
fn_image_from_mat(Probs)
grp = get(get(gcf, 'Children'), 'Children');
for im = 2:22
    grp(im).CLim = [0, alpha];
    grp(im).XLim = [Ims(1).x(1)*1e3, Ims(1).x(end)*1e3];
    grp(im).YLim = [Ims(1).z(1)*1e3, Ims(1).z(end)*1e3];
end
grp(1).Label.String = "p";
if ~strcmp(model_options.model.savepath, "")
    savefig(fullfile(model_options.model.savepath, sprintf("%s View-wise probabilities", model_options.model.savename)))
end
                
fn_image_from_mat(Fused)
grp = get(get(gcf, 'Children'), 'Children');
grp(2).CLim = [0, alpha];
grp(2).XLim = [Ims(1).x(1)*1e3, Ims(1).x(end)*1e3];
grp(2).YLim = [Ims(1).z(1)*1e3, Ims(1).z(end)*1e3];
grp(1).Label.String = "p";
if ~strcmp(model_options.model.savepath, "")
    savefig(fullfile(model_options.model.savepath, sprintf("%s Fisher fusion", model_options.model.savename)))
    saveas(gcf, fullfile(model_options.model.savepath, sprintf("%s Fisher fusion.png", model_options.model.savename)))
    save(fullfile(model_options.model.savepath, sprintf("%s Fisher fusion.mat", model_options.model.savename)), "Probs", "Fused")
end

end
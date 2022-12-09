function [Fused, Weighted, Ims] = fn_simple_mmatched_fusion(model_options)
% Performs a simple matched filter fusion. Computes the sensitivity maps
% for multiple views for a sample and TFMs from an FMC, multiplying one by
% the other pixel-by-pixel, before summing over all views to achieve a
% single fused image. Note that grain noise is assumed non-existent thus
% σ=1.

%% Compute sensitivity maps E and TFMs I
[Sens, Views] = fn_sens(model_options);

% model_options.mesh.scat.fmc_mask = 1;

% Reuse Views as Views(view).min_times were calculated alread
[Ims, ~, ~]   = fn_tfm(model_options);%, Views);

number_of_views = size(Ims, 1);
% mask = model_options.model.fusion_mask;


%% Weight TFMs using E
Weighted = Ims;
max_ = 0;
for view = 1:number_of_views
%     if isstruct(mask)
%         Weighted(view).image = Ims(view).image .* Sens(view).image;% .* mask(view).fusion_mask;
%     else
        Weighted(view).image = Ims(view).image .* Sens(view).image;% .* mask;
%     end
    if max(abs(Weighted(view).image(:))) > max_
        max_ = max(abs(Weighted(view).image(:)));
    end
end

for view = 1 : number_of_views
    Weighted(view).db_image = 20 * log10(abs(Weighted(view).image) ./ max_); 
end


%% Fused
Fused = Ims(1);
Fused.name = "Fused TFM";
Fused.image = 0;
for view = 1:number_of_views
    Fused.image = Fused.image + Weighted(view).image;
end
Fused.db_image = 20 * log10(abs(Fused.image) ./ max(abs(Fused.image(:))));

fn_image_from_mat(Fused)

if model_options.model.savepath ~= ""
    cd(model_options.model.savepath)
    filename_fig = sprintf('%s.fig', model_options.model.savename);
    filename_mat = sprintf('%s.mat', model_options.model.savename);
    savefig(filename_fig)
    save(filename_mat, 'Fused', 'Weighted')
end

end
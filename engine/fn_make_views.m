function Views = fn_make_views(Paths, unique_only)
% Creates all the unique views from paths, making sure that the same ones
% as used in arim are used, for comparison.
%
% INPUTS:
% - VIEWS : integer
%       Integer which indicates how many views are being examined. 1 means
%       direct, 2 means skip from backwall, 3 means skip from sidewall, 4
%       means skip from both backwall and sidewall.
% - Paths : array (no_paths, 1)
%       The path objects output from fn_compute_ray.
% - unique_only : logical
%       Logical test for whether to make unique views only or not. When
%       simulating, this should be 0, and when imaging this should be 1.
%
% OUTPUTS:
% - Views : struct (number_of_ims, 1)

no_paths = size(Paths, 2);

% Initialise array for storing names.
tx_rx_names = repmat("", no_paths^2, 2);
tx_rx_names_rev = repmat("", no_paths^2, 2);
% Initialise array for storing the index of the corresponding path within
% Paths array.
tx_rx_idxs = zeros(no_paths^2, 2);

ii = 1;
% Get all of the possible view names - keeping them separated by transmit
% and receive path for now - so that we can sort them into the right order.
for tx = 1:no_paths
    for rx = 1:no_paths
        tx_rx_names(ii, :) = [Paths(tx).path_info.name, Paths(rx).path_info.rev_name];
        tx_rx_idxs(ii, :) = [tx, rx];
        
        tx_rx_names_rev(ii, :) = [Paths(rx).path_info.name, Paths(tx).path_info.rev_name];

        ii = ii+1;
    end
end



% Sort view names. To be the same as arim, view names need to be sorted 
% into ascending order according to the follow hierarchy:
% 1) the total number of legs.
% 2) the maximum number of legs for transmit and receive paths.
% 3) the number of legs for the receive path.
% 4) the number of legs for the transmit path.
% 5) the lexicographic order for the transmit path.
% 6) the lexicographic order for the receive path.

% Sort criteria (6).
[~,I] = sort(tx_rx_names(:,2));
tx_rx_names(:,1) = tx_rx_names(I,1);
tx_rx_names(:,2) = tx_rx_names(I,2);
tx_rx_idxs(:,1) = tx_rx_idxs(I,1);
tx_rx_idxs(:,2) = tx_rx_idxs(I,2);

tx_rx_names_rev(:,1) = tx_rx_names_rev(I,1);
tx_rx_names_rev(:,2) = tx_rx_names_rev(I,2);

% Sort criteria (5).
[~,I] = sort(tx_rx_names(:,1));
tx_rx_names(:,1) = tx_rx_names(I,1);
tx_rx_names(:,2) = tx_rx_names(I,2);
tx_rx_idxs(:,1) = tx_rx_idxs(I,1);
tx_rx_idxs(:,2) = tx_rx_idxs(I,2);

tx_rx_names_rev(:,1) = tx_rx_names_rev(I,1);
tx_rx_names_rev(:,2) = tx_rx_names_rev(I,2);

% Sort criteria (4).
tx_legs = count(tx_rx_names, "L") + count(tx_rx_names, "T");
[~,I] = sort(tx_legs(:,1));
tx_rx_names(:,1) = tx_rx_names(I,1);
tx_rx_names(:,2) = tx_rx_names(I,2);
tx_rx_idxs(:,1) = tx_rx_idxs(I,1);
tx_rx_idxs(:,2) = tx_rx_idxs(I,2);

tx_rx_names_rev(:,1) = tx_rx_names_rev(I,1);
tx_rx_names_rev(:,2) = tx_rx_names_rev(I,2);

% Sort criteria (3).
rx_legs = count(tx_rx_names, "L") + count(tx_rx_names, "T");
[~,I] = sort(rx_legs(:,2));
tx_rx_names(:,1) = tx_rx_names(I,1);
tx_rx_names(:,2) = tx_rx_names(I,2);
tx_rx_idxs(:,1) = tx_rx_idxs(I,1);
tx_rx_idxs(:,2) = tx_rx_idxs(I,2);

tx_rx_names_rev(:,1) = tx_rx_names_rev(I,1);
tx_rx_names_rev(:,2) = tx_rx_names_rev(I,2);

% Sort criteria (2).
max_legs = max(count(tx_rx_names, "L") + count(tx_rx_names, "T"), [], 2);
[~,I] = sort(max_legs);
tx_rx_names(:,1) = tx_rx_names(I,1);
tx_rx_names(:,2) = tx_rx_names(I,2);
tx_rx_idxs(:,1) = tx_rx_idxs(I,1);
tx_rx_idxs(:,2) = tx_rx_idxs(I,2);

tx_rx_names_rev(:,1) = tx_rx_names_rev(I,1);
tx_rx_names_rev(:,2) = tx_rx_names_rev(I,2);

% Sort criteria (1).
total_legs = sum(count(tx_rx_names, "L") + count(tx_rx_names, "T"), 2);
[~,I] = sort(total_legs);
tx_rx_names(:,1) = tx_rx_names(I,1);
tx_rx_names(:,2) = tx_rx_names(I,2);
tx_rx_idxs(:,1) = tx_rx_idxs(I,1);
tx_rx_idxs(:,2) = tx_rx_idxs(I,2);

tx_rx_names_rev(:,1) = tx_rx_names_rev(I,1);
tx_rx_names_rev(:,2) = tx_rx_names_rev(I,2);



if unique_only
    tx_rx_unique_names = repmat("", no_paths^2, 2);
    tx_rx_unique_idxs = zeros(no_paths^2, 2);
    Revnames = repmat("", no_paths^2, 1);
    ii = 1;
    % Now sort out to only get the unique names.
    for view = 1:size(tx_rx_names, 1)
        forward_name = sprintf("%s - %s", tx_rx_names(view, 1), tx_rx_names(view, 2));
        if ~ismember(forward_name, Revnames)
            reverse_name = sprintf("%s - %s", tx_rx_names_rev(view, 1), tx_rx_names_rev(view, 2));
            Revnames(ii) = reverse_name;
            tx_rx_unique_names(ii, :) = tx_rx_names(view, :);
            tx_rx_unique_idxs(ii, :) = tx_rx_idxs(view, :);
            ii = ii+1;
        end
    end

    % Remove unused rows.
    tx_rx_names = tx_rx_unique_names(tx_rx_unique_names(:,1) ~= "", :);
    tx_rx_names = tx_rx_names(tx_rx_names(:,2) ~= "", :);
    tx_rx_idxs = tx_rx_unique_idxs(tx_rx_unique_idxs(:,1) ~= 0, :);
    tx_rx_idxs = tx_rx_idxs(tx_rx_idxs(:,2) ~= 0, :);
end



% Names are now sorted. We can now make the views in the right order.
no_views = size(tx_rx_idxs, 1);
Views = repmat(fn_create_view(Paths(1), Paths(1)), no_views, 1);
for view = 1:no_views
    Views(view) = fn_create_view(Paths(tx_rx_idxs(view,1)), Paths(tx_rx_idxs(view,2)));
    Views(view).name = sprintf("%s - %s", tx_rx_names(view, 1), tx_rx_names(view, 2));
end

end
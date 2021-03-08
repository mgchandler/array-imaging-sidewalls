function Views = fn_make_views(VIEWS, Paths, unique_only)
% Creates all the unique views from paths, making sure that the same ones
% as used in arim are used, for comparison. Note that arim has a sort
% algorithm to do this, for now we simply have a list of all possible views
% and reference that. Down the line, a more formal implementation should be
% considered (as in arim.ut.default_viewname_order).
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
Revnames = repmat("", no_paths^2, 1);
% Initialise array for storing the index of the corresponding path within
% Paths array.
tx_rx_idxs = zeros(no_paths^2, 2);

ii = 1;
% Get all of the possible view names - keeping them separated by transmit
% and receive path for now - so that we can sort them into the right order.
% If we only want unique names, remove duplicates.
for tx = 1:no_paths
    for rx = 1:no_paths
        % If we want unique names, check if we've already seen it.
        if unique_only
            viewname = sprintf("%s - %s", Paths(tx).path_info.name, Paths(rx).path_info.rev_name);
            if ~ismember(viewname, Revnames)
                revname = sprintf("%s - %s", Paths(rx).path_info.name, Paths(tx).path_info.rev_name);
                Revnames(ii) = revname;
                
                tx_rx_names(ii, :) = [Paths(tx).path_info.name, Paths(rx).path_info.rev_name];
                tx_rx_idxs(ii, :) = [tx, rx];
                ii = ii + 1;
            end
        % If we want all names, just get them all.
        else
            tx_rx_names(ii, :) = [Paths(tx).path_info.name, Paths(rx).path_info.rev_name];
            tx_rx_idxs(ii, :) = [tx, rx];
            
            ii = ii+1;
        end
    end
end

% Clear out unused rows.
tx_rx_names = tx_rx_names(tx_rx_names(:,1) ~= "", :);
tx_rx_names = tx_rx_names(tx_rx_names(:,2) ~= "", :);
tx_rx_idxs = tx_rx_idxs(tx_rx_idxs(:,1) ~= 0, :);
tx_rx_idxs = tx_rx_idxs(tx_rx_idxs(:,2) ~= 0, :);



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

% Sort criteria (5).
[~,I] = sort(tx_rx_names(:,1));
tx_rx_names(:,1) = tx_rx_names(I,1);
tx_rx_names(:,2) = tx_rx_names(I,2);
tx_rx_idxs(:,1) = tx_rx_idxs(I,1);
tx_rx_idxs(:,2) = tx_rx_idxs(I,2);

% Sort criteria (4).
tx_legs = count(tx_rx_names, "L") + count(tx_rx_names, "T");
[~,I] = sort(tx_legs(:,1));
tx_rx_names(:,1) = tx_rx_names(I,1);
tx_rx_names(:,2) = tx_rx_names(I,2);
tx_rx_idxs(:,1) = tx_rx_idxs(I,1);
tx_rx_idxs(:,2) = tx_rx_idxs(I,2);

% Sort criteria (3).
rx_legs = count(tx_rx_names, "L") + count(tx_rx_names, "T");
[~,I] = sort(rx_legs(:,2));
tx_rx_names(:,1) = tx_rx_names(I,1);
tx_rx_names(:,2) = tx_rx_names(I,2);
tx_rx_idxs(:,1) = tx_rx_idxs(I,1);
tx_rx_idxs(:,2) = tx_rx_idxs(I,2);

% Sort criteria (2).
max_legs = max(count(tx_rx_names, "L") + count(tx_rx_names, "T"), [], 2);
[~,I] = sort(max_legs);
tx_rx_names(:,1) = tx_rx_names(I,1);
tx_rx_names(:,2) = tx_rx_names(I,2);
tx_rx_idxs(:,1) = tx_rx_idxs(I,1);
tx_rx_idxs(:,2) = tx_rx_idxs(I,2);

% Sort criteria (1).
total_legs = sum(count(tx_rx_names, "L") + count(tx_rx_names, "T"), 2);
[~,I] = sort(total_legs);
tx_rx_names(:,1) = tx_rx_names(I,1);
tx_rx_names(:,2) = tx_rx_names(I,2);
tx_rx_idxs(:,1) = tx_rx_idxs(I,1);
tx_rx_idxs(:,2) = tx_rx_idxs(I,2);



% Names are now sorted. We can now make the views in the right order.
no_views = size(tx_rx_idxs, 1);
Views = repmat(fn_create_view(Paths(1), Paths(1)), no_views, 1);
for view = 1:no_views
    Views(view) = fn_create_view(Paths(tx_rx_idxs(tx,1)), Paths(tx_rx_idxs(rx,2)));
    Views(view).name = sprintf("%s - %s", tx_rx_names(view, 1), tx_rx_names(view, 2));
end

% if unique_only
%     
%     no_paths = size(Paths, 2);
%     View_names = repmat('L-L', no_paths^2, 1);
%     view = 1;
%     for tx = 1:no_paths
%         for rx = 1:no_paths
%             View_names(view) = sprintf("%s-%s", Paths(tx).name, Paths(rx).rev_name);
%             
%             view = view + 1;
%         end
%     end
%     
%     
%     DIRECT_IMAGING_MODES = [ ...
%                       "L-L", "L-T", "T-T" ...
%     ];
%     BACKWALL_IMAGING_MODES = [ ...
%                       "L-L", "L-T", "T-T", "LBL-L", "LBL-T", "LBT-L", "LBT-T", ...
%                       "TBL-L", "TBL-T", "TBT-L", "TBT-T", "LBL-LBL", "LBL-LBT", ...
%                       "LBL-TBL", "LBL-TBT", "LBT-LBT", "LBT-TBL", "LBT-TBT", ...
%                       "TBL-LBT", "TBL-TBT", "TBT-TBT" ...
%     ];
%     SIDEWALL_IMAGING_MODES = [ ....
%                       "L-L", "L-T", "T-T", "LSL-L", "LSL-T", "LST-L", "LST-T", ...
%                       "TSL-L", "TSL-T", "TST-L", "TST-T", "LSL-LSL", "LSL-LST", ...
%                       "LSL-TSL", "LSL-TST", "LST-LST", "LST-TSL", "LST-TST", ...
%                       "TSL-LST", "TSL-TST", "TST-TST" ....
%     ];
%     BACK_SIDEWALL_IMAGING_MODES = [....
%                       "L-L", "L-T", "T-T", "LBL-L", "LBL-T", "LBT-L", "LBT-T", ...
%                       "LSL-L", "LSL-T", "LST-L", "LST-T", "TBL-L", "TBL-T", ...
%                       "TBT-L", "TBT-T", "TSL-L", "TSL-T", "TST-L", "TST-T", ...
%                       "LBL-LBL", "LBL-LBT", "LBL-LSL", "LBL-LST", "LBL-TBL", ...
%                       "LBL-TBT", "LBL-TSL", "LBL-TST", "LBT-LBT", "LBT-LSL", ...
%                       "LBT-LST", "LBT-TBL", "LBT-TBT", "LBT-TSL", "LBT-TST", ...
%                       "LSL-LBT", "LSL-LSL", "LSL-LST", "LSL-TBT", "LSL-TSL", ...
%                       "LSL-TST", "LST-LBT", "LST-LST", "LST-TBT", "LST-TSL", ...
%                       "LST-TST", "TBL-LBT", "TBL-LST", "TBL-TBT", "TBL-TST", ...
%                       "TBT-LST", "TBT-TBT", "TBT-TST", "TSL-LST", "TSL-TST", ...
%                       "TST-TST" ...
%     ];
% 
% 
% 
%     if VIEWS == 1
%         modes = DIRECT_IMAGING_MODES;
%         Number_of_ims = 3;
%     elseif VIEWS == 2
%         modes = BACKWALL_IMAGING_MODES;
%         Number_of_ims = 21;
%     elseif VIEWS == 3
%         modes = SIDEWALL_IMAGING_MODES;
%         Number_of_ims = 21;
%     elseif VIEWS == 4
%         modes = BACK_SIDEWALL_IMAGING_MODES;
%         Number_of_ims = 55;
%     else
%         error('fn_make_views: Invalid VIEWS parameter')
%     end
% 
%     ii = 1;
% 
% 
% 
%     Views = repmat(fn_create_view(Paths(1), Paths(1)), Number_of_ims, 1);
% 
%     for tx = 1 : size(Paths, 2)
%         for rx = 1 : size(Paths, 2)
%             thisPathName = Paths(tx).path_info.name;
%             thisRevPathName = reverse(Paths(rx).path_info.name);
%             name = strcat(thisPathName, "-", thisRevPathName);
%             if any(modes == name)
%                 Views(ii) = fn_create_view(Paths(tx), Paths(rx));
%                 Views(ii).name = name;
%                 ii = ii+1;
%             end
%         end
%     end
%     
% % Make all views
% else
%     ii = 1;
% 
% 
% 
%     Views = repmat(fn_create_view(Paths(1), Paths(1)), size(Paths, 1)^2, 1);
% 
%     for tx = 1 : size(Paths, 2)
%         for rx = 1 : size(Paths, 2)
%             thisPathName = Paths(tx).path_info.name;
%             thisRevPathName = reverse(Paths(rx).path_info.name);
%             name = strcat(thisPathName, "-", thisRevPathName);
%             Views(ii) = fn_create_view(Paths(tx), Paths(rx));
%             Views(ii).name = name;
%             ii = ii+1;
%         end
%     end
%     
% end

end
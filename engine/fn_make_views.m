function Views = fn_make_views(VIEWS, Paths, varargin)
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
%
% OUTPUTS:
% - Views : struct (number_of_ims, 1)

DIRECT_IMAGING_MODES = [ ...
                  "L-L", "L-T", "T-T" ...
];
BACKWALL_IMAGING_MODES = [ ...
                  "L-L", "L-T", "T-T", "LBL-L", "LBL-T", "LBT-L", "LBT-T", ...
                  "TBL-L", "TBL-T", "TBT-L", "TBT-T", "LBL-LBL", "LBL-LBT", ...
                  "LBL-TBL", "LBL-TBT", "LBT-LBT", "LBT-TBL", "LBT-TBT", ...
                  "TBL-LBT", "TBL-TBT", "TBT-TBT" ...
];
SIDEWALL_IMAGING_MODES = [ ....
                  "L-L", "L-T", "T-T", "LSL-L", "LSL-T", "LST-L", "LST-T", ...
                  "TSL-L", "TSL-T", "TST-L", "TST-T", "LSL-LSL", "LSL-LST", ...
                  "LSL-TSL", "LSL-TST", "LST-LST", "LST-TSL", "LST-TST", ...
                  "TSL-LST", "TSL-TST", "TST-TST" ....
];
BACK_SIDEWALL_IMAGING_MODES = [....
                  "L-L", "L-T", "T-T", "LBL-L", "LBL-T", "LBT-L", "LBT-T", ...
                  "LSL-L", "LSL-T", "LST-L", "LST-T", "TBL-L", "TBL-T", ...
                  "TBT-L", "TBT-T", "TSL-L", "TSL-T", "TST-L", "TST-T", ...
                  "LBL-LBL", "LBL-LBT", "LBL-LSL", "LBL-LST", "LBL-TBL", ...
                  "LBL-TBT", "LBL-TSL", "LBL-TST", "LBT-LBT", "LBT-LSL", ...
                  "LBT-LST", "LBT-TBL", "LBT-TBT", "LBT-TSL", "LBT-TST", ...
                  "LSL-LBT", "LSL-LSL", "LSL-LST", "LSL-TBT", "LSL-TSL", ...
                  "LSL-TST", "LST-LBT", "LST-LST", "LST-TBT", "LST-TSL", ...
                  "LST-TST", "TBL-LBT", "TBL-LST", "TBL-TBT", "TBL-TST", ...
                  "TBT-LST", "TBT-TBT", "TBT-TST", "TSL-LST", "TSL-TST", ...
                  "TST-TST" ...
];



if VIEWS == 1
    modes = DIRECT_IMAGING_MODES;
    Number_of_ims = 3;
elseif VIEWS == 2
    modes = BACKWALL_IMAGING_MODES;
    Number_of_ims = 21;
elseif VIEWS == 3
    modes = SIDEWALL_IMAGING_MODES;
    Number_of_ims = 21;
elseif VIEWS == 4
    modes = BACK_SIDEWALL_IMAGING_MODES;
    Number_of_ims = 55;
else
    error('fn_make_views: Invalid VIEWS parameter')
end

i = 1;

% If we need to use scat_info while creating views, collect this.
if nargin > 2
    scat_info = varargin{1};
    Views = repmat(fn_create_view(Paths(1), Paths(1), scat_info), Number_of_ims, 1);

    for tx = 1 : size(Paths, 2)
        for rx = 1 : size(Paths, 2)
            thisPathName = Paths(tx).path_info.name;
            thisRevPathName = reverse(Paths(rx).path_info.name);
            name = strcat(thisPathName, "-", thisRevPathName);
            if any(modes == name)
                Views(i) = fn_create_view(Paths(tx), Paths(rx), scat_info);
                Views(i).name = name;
                i = i+1;
            end
        end
    end
% If we do not need scat_info, then do not pass this through.
else
    Views = repmat(fn_create_view(Paths(1), Paths(1)), Number_of_ims, 1);
    for tx = 1 : size(Paths, 2)
        for rx = 1 : size(Paths, 2)
            thisPathName = Paths(tx).path_info.name;
            thisRevPathName = reverse(Paths(rx).path_info.name);
            name = strcat(thisPathName, "-", thisRevPathName);
            if any(modes == name)
                Views(i) = fn_create_view(Paths(tx), Paths(rx));
                Views(i).name = name;
                i = i+1;
            end
        end
    end
end

end
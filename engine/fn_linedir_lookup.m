function dir_out = fn_linedir_lookup(theta, mode, rad, npw)
% Simple lookup of directivity measured in an FE simulation when centre
% frequency was 5MHz and element width was 0.0mm operating in contact with
% an aluminium block. Inputs may be vectors but are required to have same
% length, except for mode, which should always be a single integer. Returns
% directivities of equal length to inputs.

assert(length(mode) == 1, 'fn_linedir_lookup: mode should have length 1.')

% Read table
filename = sprintf('D_rdm1_%dnpw.dir_32interp.lookup.dat', npw);
dir = readtable(filename);

% Get unique radii for which we have directivity values.
rmags_read = dir.rmag([true; (diff(dir.rmag)~=0)]);
numPhi = height(dir) / length(rmags_read);

% Generate rmags and angles which have uniform spacing
rmags = linspace(min(rmags_read), max(rmags_read), length(rmags_read));
phi = linspace(min(table2array(dir(dir.rmag==rmags(1), 2))), max(table2array(dir(dir.rmag==rmags(1), 2))), numPhi);

% Check that all radii have the same number of angles, and that the
% uniformly spaced angles are similar enough to the actual angles to use
% for the rest of the analysis.
assert(all(abs(rmags_read' - rmags) < 1e-10), 'fn_linedir_lookup: rmags not similar enough.')
for r = 1:length(rmags)
    numTheta = size(table2array(dir(dir.rmag == rmags_read(r), 2)), 1);
    assert(numTheta == numPhi, 'fn_linedir_lookup: incorrect number of phi in directivity table')
    assert(all(abs(table2array(dir(dir.rmag == rmags_read(r), 2))' - phi) < 1e-10), sprintf('fn_linedir_lookup: phi not similar enough for rmag %.2f', rmags(r)))
end

% We have the same number of phi for each radial distance. Convert the
% directivity values to a matrix to interpolate.
dL = reshape(table2array(dir(:, 3)), numPhi, length(rmags));
dS = reshape(table2array(dir(:, 4)), numPhi, length(rmags));

if ~isa(dL, 'double')
    dL = str2double(strrep(strrep(dL, '(', ''), ')', ''));
end
if ~isa(dS, 'double')
    dS = str2double(strrep(strrep(dS, '(', ''), ')', ''));
end

% Interpolate within directivity data.
if mode
    dir_out = interp2(rmags, phi, dS, rad, theta, 'cubic', 0);
else
    dir_out = interp2(rmags, phi, dL, rad, theta, 'cubic', 0);
end

end
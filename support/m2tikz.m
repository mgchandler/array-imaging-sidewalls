function m2tikz(filename, varargin)
% Wrapper for matlab2tikz which solves the problems I have had. Leaves all
% text in the file, sets the figure size to 80% of linewidth, and includes
% the option for a relative path to figure location.

if ~isempty(findobj('type', 'figure'))
    filename = char(filename);
    cleanfigure('pruneText', false)
    if nargin == 1
        matlab2tikz(filename, 'width', '0.8\linewidth', 'mathmode', false)
    else
        path = char(varargin{1});
        matlab2tikz(filename, 'width', '0.8\linewidth', 'mathmode', false, 'relativeDataPath', path)
    end
else
    return
end
end
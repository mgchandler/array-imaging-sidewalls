function m2tikz(filename, varargin)
% Wrapper for matlab2tikz which solves the problems I have had. Leaves all
% text in the file, sets the figure size to 80% of linewidth, and includes
% the option for a relative path to figure location.
%
% Syntax:
%   m2tikz(filename)
%   m2tikz(filename, path)
%   m2tikz(___, Name, Value)
%
% Parameters:
%   filename : char or string
%   path: char or string
%   Name-Value Arguments:
%       Passed directly into matlab2tikz - see their documentation.

if ~isempty(findobj('type', 'figure'))
    filename = char(filename);
    cleanfigure('pruneText', false)
    % Do the default process
    if nargin == 1
        matlab2tikz(filename, 'width', '0.8\linewidth', 'mathmode', false)
    % More options provided on top, unpack using defaults.
    elseif nargin > 1
        % Make a structure which contains the default arguments. Then
        % overwrite them later if they are specified.
        args = struct('width', '0.8\linewidth', 'mathmode', false);
        % Path may or may not be listed with its argument handle name -
        % check if it is or not.
        if mod(nargin, 2) == 0
            % Path must be listed without argument name
            args.relativeDataPath = char(varargin{1});
            start = 1;
        else
            % Else path must be listed with argument name
            start = 0;
        end
        for arg = 1:(nargin-start-1)/2
            args.(varargin{1 + start + 2*(arg-1)}) = char(varargin{1 + start + 2*(arg-1) + 1});
        end
        unpack = {};
        names = fieldnames(args);
        vals = struct2cell(args);
        for arg = 1:size(names, 1)
            unpack = {unpack{:}, names{arg}, vals{arg}};
        end
        matlab2tikz(filename, unpack{:})
    end
else
    return
end
end
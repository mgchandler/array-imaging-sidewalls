function h = fn_acl(x, y, im, varargin)

if numel(x)*numel(y) ~= numel(im)
    error("fn_acf: x and y should be broadcastable with im")
end
if nargin > 3
    threshold = varargin{1};
else
    threshold = 1/exp(1);
end

dx = x(2) - x(1);
dy = y(2) - y(1);

q = xcorr2(im);
q = q ./ max(q(:));

X = [1:size(q, 1)]*dx;
Y = [1:size(q, 2)]*dy;
X = X - mean(X);
Y = Y - mean(Y);

gr_e = q >= threshold;
q(~gr_e) = 0;

[meshX, meshY] = meshgrid(X, Y);
R = sqrt(meshX.^2 + meshY.^2);
h = max(R(gr_e));

end
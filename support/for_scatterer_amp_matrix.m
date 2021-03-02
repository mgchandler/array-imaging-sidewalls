sdh = fn_scat_info( ...
    "sdh", ...
    1.0e-3, ...
    6320.0/5.0e+6, ...
    3130.0/5.0e+6, ...
    0, ...
    [[0, 0, 0]], ...
    'ang_pts_over_2pi', 60 ...
);

max_ = 0;
for mode = 1:4
    matrix = sdh.matrix(mode);
    if max(abs(matrix.matr(:))) > max_
        max_ = max(abs(matrix.matr(:)));
    end
end

LL = abs(sdh.matrix(1).matr) / max_;
LT = abs(sdh.matrix(2).matr) / max_;
TL = abs(sdh.matrix(3).matr) / max_;
TT = abs(sdh.matrix(4).matr) / max_;

theta = linspace(-pi, pi, 120);

subplot(2,2,1)
imagesc(theta, theta, LL)
caxis([0, 1])

subplot(2,2,2)
imagesc(theta, theta, LT)
caxis([0, 1])

subplot(2,2,3)
imagesc(theta, theta, TL)
caxis([0, 1])

subplot(2,2,4)
imagesc(theta, theta, TT)
caxis([0, 1])

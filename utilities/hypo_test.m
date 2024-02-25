%% hypothesis testing for beta(s): C * beta(s) = c(s)
% Input:
%       C: full rank matrix in the hypothesis, (r,p)
%       c: each row is a vecotr with values function c(s) at m grids, (r * m)

% Output:
function [pval, Tn0, Tn_R] = hypo_test(x, ally, t, tau, hx, hy, h, ...
    betaest, gest, dgest, h2, h3, C, c, R, smooth)
if ~exist('smooth', 'var')
    smooth = false;
end

%% Obtain the test statistics Tn0 using original dataset
[n,p] = size(x);
ystar = ally - gest;
m = size(ystar, 2);
etaest = zeros(n, m);
for i=1:n
    for s=1:m
        % solve eq(14), view Wms=1 to estimate eta_i(s)
        etaest(i, s) = locallinear0(1, t(s), h2, t', ystar(i, :)');
    end
end
Tn0 = global_Tn(C, c, betaest, dgest, x, etaest, t, tau, hx, h3);

%% Get the estimators under H0
idx = sum(C, 1) == 1;
x_H0 = x; x_H0(:, idx) = []; % (n, p0)
beta0 = betaest; beta0(idx, :) = []; % (p0, m)
[beta_H0, gest_H0, dgest_H0, h1] = get_estimate(x_H0, ally, t, tau, ...
    hx, hy, h, beta0, gest, dgest, smooth);
betaest_H0 = betaest; betaest_H0(~idx, :) = beta_H0;

%% Bootstrap & get the test statistics Tn_R under H0
rng(2021);
all_zeta = randn(n,R);
Tn_R = zeros(1, R);
res = ally - gest_H0;
for r = 1:R
    zeta = all_zeta(:, r);
    yboot = gest_H0 + zeta .* res;
    [betaest_boot, gest_boot, dgest_boot, ~] = get_estimate(x, yboot, t, ...
        tau, hx, hy, h, betaest_H0, gest_H0, dgest_H0, smooth, h1);
    ystar_boot = yboot - gest_boot;
    etaest_boot = zeros(n, m);
    for i=1:n
        for s=1:m
            etaest_boot(i, s) = locallinear0(1, t(s), h2, t', ystar_boot(i, :)');
        end
    end
    Tn_R(r) = global_Tn(C, c, betaest_boot, dgest_boot, x, etaest_boot, t, tau, hx, h3);
end

%% Calculate p-value using chi2 approximation
k1 = mean(Tn_R);
k2 = var(Tn_R);
k3 = mean((Tn_R - k1).^3);
a1 = k3 / ( 4 * k2);
a2 = 8 * k2^3 / (k3^2);
a3 = k1 - 2 * k2^2 / k3;
pval = chi2cdf((Tn0 - a3) / a1, a2, 'upper');

end

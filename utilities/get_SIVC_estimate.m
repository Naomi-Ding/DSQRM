%% Estimation for SIVC
function [betaest, gest] = get_SIVC_estimate(x, ally, t, hx, hy, h, beta0)

[n,p] = size(x);
m = size(ally, 2);

%% Step 1. Estimate Beta %%
% disp("Estimate beta(s)");
betaest = zeros(p,m);
for s=1:m
    betaest(:, s) = getBeta(x, ally, m, t, s, hx, hy, h, beta0(:, s));
end

%% Step 2. Estimate g %%
% disp("Estimate g")
gest = zeros(n,m);
h2 = 0.014;
allxb = x * betaest; % X.T * beta(s), (n,m)
for i=1:n
    for s=1:m
        xb0 = x(i, :)*betaest(:, s);
        gest(i, s) = locallinear0(1, xb0, h2, allxb(:), ally(:));
    end
end

end
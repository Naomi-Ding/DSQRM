%% Inference for g(z): SCB for g(z)

function [Cb_g] = inference_g(x, ally, betaest, gest, dgest, tau, h1, R, a)

rng(2021);
[n,m] = size(ally);
all_zeta = randn(n,m,R);
Cg_r = zeros(1,R);
for r = 1:R
    zeta = all_zeta(:,:,r);
    Cg_r(r) = SCB_g(x, ally, betaest, gest, dgest, tau, h1, zeta);
end

tmp = sort(abs(Cg_r));
Cb_g = tmp(ceil(R*(1 - a)));

end

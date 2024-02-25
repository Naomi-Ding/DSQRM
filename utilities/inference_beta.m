%% Inference for beta(s): SCB for beta_l(s), l = 1,..,p; SCR for beta(s)
% Input:
%       x: covariates of interest, (n,p)
%       ystar: residuals after deducting g from the response, (n, m)
%       tau: the targeted quantile level, scalar
%       R: number of bootstraps in estimating SCB & SCR for beta(s) & g
%       a: alpha, significance level used in the inference
%       varargin: a cell array containing (h2, h3),
%                 if not provided, use the default settings
%               bandwidth: (h2, h3)
% Output:
%       SCB_beta: simultaneous confidence bands for beta_l(s), l=1,...,p, (p,1)
%       SCR_beta: simultaneous confidence region for beta(s), scalar
function [SCB_beta, SCR_beta, h2, h3] = inference_beta(x, ystar, betaest, dgest, t, tau, hx, R, a, varargin)
rng(2021);

[n,p] = size(x);
m = size(ystar, 2);

if nargin > 8
    h2 = varargin{1};
else
    h2 = cvh2(t, ystar);
end
etaest = zeros(n, m);
for i=1:n
    for s=1:m
        % solve eq(14), view Wms=1 to estimate eta_i(s)
        etaest(i, s) = locallinear0(1, t(s), h2, t', ystar(i, :)');
    end
end

if nargin > 9
    h3 = varargin{2};
else
    f_SCV = @(h3) SCV(etaest, h3);
    h3 = fminbnd(f_SCV, 0, 1);
end

C = ones(1,p); c = zeros(1,m);
[~, feta_0, all_A, all_inv_meanB] = global_Tn(C, c, betaest, dgest, x, etaest, t, tau, hx, h3);
Cp_b = zeros(p,R);
Cr_b = zeros(R,1);
all_zeta = randn(n,R);
for r = 1:R
    zeta = all_zeta(:,r);
    G = zeros(p,m);
    X = zeros(m,1);
    for s = 1:m
        A = all_A(:,:,s); % (n,p), Ai'
        inv_meanB = all_inv_meanB(:,:,s); % (p,p), inverse of mean of Bi
        sumA = (zeta.* A)' * phi(etaest(:,s), tau); % (p,1), sum(zeta_i*Ai*phi_etai)
        meanA = sumA ./ sqrt(n); % (p,1)
        G(:,s) = feta_0(s) .* inv_meanB * meanA;
        X(s) = feta_0(s)^2 .* meanA' * (inv_meanB' * inv_meanB) * meanA; % scalar
    end
    Cp_b(:,r) = max(G,[],2); % (p,1), SCB for beta_l(s)
    Cr_b(r) = max(X); % scalar, SCR for beta(s) based on MELE
end

% SCB based on MELEs
sort_Cp_b = sort(Cp_b,2);
CB_b = sort_Cp_b(:,ceil(R*(1-a)));
SCB_beta = CB_b / sqrt(n);

% SCR
tmp = sort(Cr_b);
SCR_beta = tmp(ceil(R*(1-a)));

end
%% Main function for DSQRM: input image, output estimated parameters
% Input:
%       x: covariates of interest, (n,p)
%       v: pixel intensities of images, (n, 1) cell array, each cell contains N_i pixels
%       m: number of points for estimated densities of the images
%       tau_set: the targeted quantile levels, (1, ntau)
%       varargin: a cell array containing (hx, hy, h, beta0, g0, dg0, R, a, [h2, h3, C ,c]),
%                 if not provided, use the default settings
%               bandwidth: (hx, hy, h)
%               initial values of (beta(s), g(x^T beta(s))), dg(x^T beta(s))
%               R: number of bootstraps in estimating SCB & SCR for beta(s) & g
%               a: alpha, significance level used in the inference
%               bandwidth: (h2, h3), used in the inference procedure
%               C & c: the hypothesis to be tested on beta(s),
%                       if not provided, test the significance of the last covariate
% Output:
%       fhat: estimated denstiies, (n, N)
%       f_support: support of estimated densities, (n, N)
%       hf: bandwidth for density estimators, (1, n)
%       ally: LQD representation of fhat, (n, m)
%       all_betaest: estimated functional coefficients at targeted quantile levels, (p,m,ntau)
%       all_gest: estimated link function at target quantile levels, (n,m,ntau)
%       all_gest_inv: inversed transformation of estimated gest, (n,m,ntau)
%       **Optional**:
%       all_Cb_beta: simultaneous confidence bands for estimated coefficients beta_l(s), l=1,...,p; (p,ntau)
%       all_Cr_beta: simultaneous confidence region for estimated coefficient functions beta(s), (1, ntau)
%       all_Cb_g: simultaneous confidence band for estimated link function gest(\cdot), (1, ntau)
%       all_pvals: p-values of the hypothesis testing procedures at the targeted quantile levels, (1, ntau)

function [fhat, f_support, hf, ally, all_betaest, all_gest, all_dgest, all_gest_inv, ...
    all_Cb_beta, all_Cr_beta, all_Cb_g, all_pvals] = DSQRM(x, v, m, tau_set, ...
    hx, hy, h, h1, beta0, g0, dg0, smooth, verbose, SCB, Hypo_idx, R, a, h2, h3)
% addpath('./utilities')
[n,p] = size(x);
if ~exist('verbose', 'var')
    verbose = 0;
end
tic;
%% Step 1. Density Estimation & Obtain LQD
disp("Step 1. Extract density estimators & LQD representations");
[ally, t, hf, fhat, f_support] = get_lqd_representation(v, m);
toc;

rho = std(x(:));
if ~exist('hx', "var") || isempty(hx)
    hx = n^(-1/3)*rho^0.5;   %cn^(-1/3)
end

if ~exist('hy', "var") || isempty(hy)
    hy = n^(-1/5)*rho^0.5;   %cn^(-1/5)
end

if ~exist('h', "var") || isempty(h)
    h = m^(-1/5)*0.3;
end

if ~exist('beta0', "var") || isempty(beta0)
    if ~exist('h1', "var") || isempty(h1)
        [beta0, g0, dg0] = get_initials(x, ally, t, hx, hy, h, 0.15);
    else
        [beta0, g0, dg0] = get_initials(x, ally, t, hx, hy, h, h1);
    end
elseif ~exist('g0', "var") || ~exist('dg0', "var") ||isempty(g0) || isempty(dg0)
    %% initials for g(xb) & dg(xb)
    g0 = zeros(n,m);
    dg0 = zeros(n,m);
    allxb0 = x * beta0;
    if ~exist('h1', "var") || isempty(h1)
        for i = 1:n
            for s = 1:m
                xb0 = x(i,:) * beta0(:,s);
                [g0(i,s), dg0(i,s)] = locallinear1(1, xb0, 0.15, allxb0(:), ally(:));
            end
        end
    else
        for i = 1:n
            for s = 1:m
                xb0 = x(i,:) * beta0(:,s);
                [g0(i,s), dg0(i,s)] = locallinear1(1, xb0, h1, allxb0(:), ally(:));
            end
        end
    end
end

if ~exist('smooth', "var")
    smooth = false;
end


%% Step 2. Get the estimators of our method
disp("Step 2. Estimate the functional coefficients & link function");
ntau = length(tau_set);
all_betaest = zeros(p,m,ntau);
all_gest = zeros(n,m,ntau);
all_dgest = zeros(n,m,ntau);
all_gest_inv = zeros(n, m, ntau);
all_h1 = zeros(1, ntau);
for idx = 1:ntau
    tau = tau_set(idx);
    disp(compose("tau = %.1f\n", tau));
    if ~exist('h1', "var") || isempty(h1)
        [betaest, gest, dgest, h1] = get_estimate(x, ally, t, tau, hx, hy, h, beta0, g0, dg0, smooth, verbose);
    else
        [betaest, gest, dgest, h1] = get_estimate(x, ally, t, tau, hx, hy, h, beta0, g0, dg0, smooth, verbose, h1);
    end
    % Get the inverse of estimated link
    gest_inv = MakeDENsample(gest);
    all_betaest(:,:,idx) = betaest;
    all_gest(:,:,idx) = gest;
    all_dgest(:,:,idx) = dgest;
    all_h1(idx) = h1;
    all_gest_inv(:,:,idx) = gest_inv;
    toc;
end

if ~exist('SCB', "var")
    SCB = false;
end

if ~exist('Hypo_idx', 'var')
    Hypo_idx = 0;
end

if SCB || Hypo_idx
    if ~exist('R', "var")
        R = 200;
    end
    if ~exist('a', 'var')
        a = 0.05;
    end
    if ~exist('h3', 'var') || isempty(h3)
        h3 = 0.6455;
    end
    
    if SCB
        %% Step 3. Confidence Bands for the functional coefficients & link function
        disp("Step 3. Confidence Bands for the functional coefficients & link function");
        all_Cb_beta = zeros(p, ntau);
        all_Cr_beta = zeros(1, ntau);
        all_Cb_g = zeros(1, ntau);
        for idx = 1:ntau
            disp(compose("SCB for tau = %.1f\n", tau));
            tau = tau_set(idx);
            ystar = ally - all_gest(:,:,idx);
            if ~exist('h2', "var") || isempty(h2)
                h2 = cvh2(t, ystar);
            end
            [SCB_beta, SCR_beta, h2, h3] = inference_beta(x, ystar, all_betaest(:,:,idx), ...
                all_dgest(:,:,idx), t, tau, hx, R, a, h2, h3);
            
            all_Cb_beta(:, idx) = SCB_beta;
            all_Cr_beta(idx) = SCR_beta;
            Cb_g = inference_g(x, ally, all_betaest(:,:,idx), all_gest(:,:,idx),...
                all_dgest(:,:,idx), tau, all_h1(idx), R, a);
            all_Cb_g(idx) = Cb_g;
            toc;
        end
    end
    
    if Hypo_idx
        %% Step 4. Hypothesis Testing on beta(s)
        disp("Step 4. Hypothesis Testing on beta(s)");
        C = zeros(1,m); C(Hypo_idx) = 1; % test the significance of the covariate with idx
        r = size(C, 1);
        if ~exist('c', "var") || isempty(c)
            c = zeros(r, m);
        end
        
        all_pvals = zeros(1, ntau);
        for idx = 1:ntau
            tau = ta_set(idx);
            ystar = ally - all_gest(:,:,idx);
            if ~exist('h2', "var") || isempty(h2)
                h2 = cvh2(t, ystar);
            end
            [pval, ~, ~] = hypo_test(x, ally, t, tau, hx, hy, h, ...
                all_betaest(:,:,idx), all_gest(:,:,idx), all_dgest(:,:,idx), ...
                h2, h3, C, c, R);
            all_pvals(idx) = pval;
        end
        toc;
    end
    
end


end

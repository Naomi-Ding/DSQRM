%% figure 3: estimators & SCB on simulated dataset
%%% settings
clear all;
addpath('./utilities')

n = 200;
p = 2;
m = 100;
N = 1000;
nsimu = 200; % only generate 200 datasets to replicate the n-th set
tau = 0.5;
nth = 159;
smooth = 1;
SCB = 1;
Hypo_idx = 0;
R = 500;
a = 0.05;
h1 = 0.014;
h2 = 0.02;
h3 = 1;

%% (1) Data generation: sampling from quantile functions with AR(1) correlation
tic;
disp("Generate Simulated Dataset");
ratio = 1;
rho = 0.6;
lambda = [0.01 0.05];
hx = n^(-1/3)*rho^0.5;   %cn^(-1/3)
hy = n^(-1/5)*rho^0.5;   %cn^(-1/5)
h = m^(-1/5)*0.3;
muX = zeros(1, p);
SigmaX = zeros(p, p);
for i=1:p
    for j=1:p
        SigmaX(i, j) = rho^(abs(i-j));
    end
end
t = linspace(0,1,m);
betaf = @(s) [1+s.^2; ratio * (1-s).^2]; % p x m
beta0 = @(s) betaf(s) ./ vecnorm(betaf(s)); % p x m
beta0 = beta0(t);
g = @(u) sin(2*pi*u).*(u+0.5);
dg = @(u) 2*pi*cos(2*pi*u).*(u+1/2) + sin(2*pi*u);

% Step 1. Data Generation
rng(2021);
all_eps = randn(n,2,nsimu);
% generate uniformly distributed AR(1) data
k = 2;
all_ei = binornd(1, 1/2, n, N, nsimu) / k;
all_U = zeros(n, N+1, nsimu);
for j = 2:(N+1)
    all_U(:,j,:) = all_U(:,j-1,:)/k + all_ei(:,j-1,:);
end
all_U = all_U(:,2:end,:);
for nn = 1:nsimu
    % fprintf("nn = %d\n", nn);
    x0 =  mvnrnd(muX, SigmaX, n);
    if nn == nth
        x = normcdf(x0);
        % simu_x(:, :, nn) = x;
        U = all_U(:,:,nn); % n x mi
        eps = all_eps(:,:,nn); % n x 2
        v = Q_error(U, ratio, 0.5, x, lambda, eps);
        % all_v(:,:,nn) = v;
        break
    end
end
toc;


%% (2) Density Estimation & Obtain LQD
disp("Step 1. Extract density estimators & LQD representations");
[ally, t, hf, fhat, f_support] = get_lqd_representation(v, m);
toc;

%% (3) Obtain the estimators
disp("Step 2. Estimate the functional coefficients & link function");
g0 = g(x*beta0);
dg0 = dg(x*beta0);
% disp(compose("tau = %.1f\n", tau));
[betaest, gest, dgest, h1] = get_estimate(x, ally, t, tau, hx, hy, h, beta0, g0, dg0, smooth, 1, h1);
% Get the inverse of estimated link
% gest_inv = MakeDENsample(gest);
toc;

%% (4) SCB
disp("Step 3. Confidence Bands for the functional coefficients & link function");
ystar = ally - gest;
[SCB_beta, ~, h2, h3] = inference_beta(x, ystar, betaest, ...
    dgest, t, tau, hx, R, a, h2, h3);
% Cb_g = inference_g(x, ally, betaest, gest, dgest, tau, h1, R, a);
Cb_g = 0.1871;
toc;

%% (5) plot
betaL = betaest - SCB_beta;
betaU = betaest + SCB_beta;
gL = gest - Cb_g;
gU = gest + Cb_g;

figure(3); clf;
for k = 1:p
    subplot(1,4,k); plot(t, betaest(k,:), 'k:', 'LineWidth',1.5); hold on;
    plot(t, betaL(k,:), 'k--', t, betaU(k,:), 'k--','LineWidth',1.5);
    plot(t, beta0(k,:), 'k-', 'LineWidth',1.5);
    %     legend('estimation', 'lower bound', 'upper bound','truth');
    %     title(sprintf('beta_%d(s)', k));
    xlabel('s', 'Fontsize',20);
end
subplot(1,4,1); title('\beta_1(s)', 'Fontsize',20);
subplot(1,4,2); title('\beta_2(s)', 'Fontsize',20);

allxb = x * betaest;
g0 = g(allxb);
xb_plot = allxb(:); [xb_plot, index] = sort(xb_plot);
g0_plot = g0(:); g0_plot = g0_plot(index);
% g_plot = gest(:);
g_plot = gest(:);
g_plot = g_plot(index);
gL_plot = gL(:); gL_plot = gL_plot(index);
gU_plot = gU(:); gU_plot = gU_plot(index);
% figure;
subplot(1,2,2);
plot(xb_plot,  g_plot, 'k:',xb_plot, gL_plot,'k--',xb_plot,gU_plot,'k--', 'LineWidth',1.5);
hold on; plot(xb_plot, g0_plot, 'k-', 'LineWidth',1.5);
% legend('estimation','lower bound','upper bound','truth');
% title('Estimated g(\cdot)');
title('g(x^T\beta(s))', 'Fontsize',20);xlabel('x^T\beta(s)', 'Fontsize',20);

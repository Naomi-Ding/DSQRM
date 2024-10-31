function [simu_x, simu_g, all_v, beta, g, dg, t, hx, hy, h] = simu_data_gen(n, p, m, N, nsimu, ratio)
%% Simulation Settings
% ratio = 1;
% n = 100;
% p = 2;
% m = 100;
% N = 500;
% nsimu = 200; % number of simulation
% % tau = 0.5;
if nargin < 6
    ratio = 1;
end

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

beta0 = @(s) [1+s.^2; ratio * (1-s).^2]; % p x m
beta = @(s) beta0(s) ./ vecnorm(beta0(s)); % p x m
g = @(u) sin(2*pi*u).*(u+0.5);
dg = @(u) 2*pi*cos(2*pi*u).*(u+1/2) + sin(2*pi*u);
t = linspace(0,1,m);

%% Step 1. Data Generation
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

all_v = zeros(n,N,nsimu);
simu_x = zeros(n,p,nsimu);
for nn = 1:nsimu
    % fprintf("nn = %d\n", nn);
    x0 =  mvnrnd(muX, SigmaX, n);
    x = normcdf(x0);
    simu_x(:, :, nn) = x;
    U = all_U(:,:,nn); % n x mi
    eps = all_eps(:,:,nn); % n x 2
    v = Q_error(U, ratio, 0.5, x, lambda, eps);
    all_v(:,:,nn) = v;
end


simu_g = zeros(n,m,nsimu);
for nn = 1:nsimu
    x = simu_x(:,:,nn); % n x p
    xb = x * beta(t); % n x m
    simu_g(:,:,nn) = g(xb);
end

end
%% Estimation for Distribution-on-scalar-Single-index-Quantile-Regression-Model

function [betaest, gest, dgest, h1] = get_estimate(x, ally, t, tau, ...
    hx, hy, h, beta0, g0, dg0, smooth, varargin)
[n,p] = size(x);
m = size(ally, 2);

if ~exist('smooth', 'var')
    smooth = false;
end

%% Step 1. Estimate Beta %%
% warning('off','all')
options = optimset('MaxFunEvals',100, 'MaxIter',100, 'Display', 'none');
g_method = 'fminsearch';
% disp("Estimate beta(s)")
betaest = zeros(p,m);
for s= 1:m
    % if mod(s, 10) == 0
    %     disp(compose("s = %d\n", s));
    % end
    G0 = [g0(:,s), dg0(:,s)]; % (n,2), starting point of (g, dg)
    funBeta = @(beta) ELLR(x, ally, beta, m, t, s, hx, hy, h, tau, G0, g_method);
    
    % --- use fminsearch, then normalize --- %
    tmp = fminsearch(funBeta, beta0(:,s), options); % (p,1)
    tmp = tmp/sqrt(sum(tmp.^2)); % norm=1
    betaest(:,s) = tmp;
end

% smooth the curves
if smooth
    betaest_smooth = zeros(p,m);
    for s = 1:m
        betaest_smooth(:,s) = locallinear0(p, t(s), 0.05, t', betaest');
    end
    betaest_smooth = betaest_smooth ./ sqrt(sum(betaest_smooth.^2,1));
    betaest = betaest_smooth;
end

%% Step 2. Estimate g %%
% disp("Estimate g")
gest = zeros(n,m);
dgest = zeros(n,m);
% (1) Compute h1 by CV
allxb = x * betaest; % X.T * beta(s), (n,m)
if ~isempty(varargin)
    h1 = varargin{1};
else
    h1 = cvh1(x, betaest, allxb, ally);
end
% (2) Estimate g(X.T*Beta(s))
for i = 1:n % for object i
    for s = 1:m % at grid s
        xb0 = x(i,:) * betaest(:,s); % xi.T * beta(sm), scalar
        xb = allxb(:); % (n*m,1) vector
        G0 = [g0(:,s), dg0(:,s)];
        %         options = optimset('MaxFunEvals',10000, 'MaxIter',10000, 'Display', 'off');
        funG = @(G) minG(G, xb, xb0, ally, m, h1, tau);
        G_est = fminsearch(funG, G0(i,:), options);
        gest(i,s) = G_est(1);
        dgest(i,s) = G_est(2);
    end
end

if smooth
    % smooth g
    gest_smooth = zeros(n,m);
    hg = 0.05;
    hdg = hg;
    xb = x * betaest;
    for i = 1:n
        %     i
        for s = 1:m
            xb0 = xb(i,s); % scalar
            tmp1 = zeros(n,1);
            tmp2 = zeros(n,1);
            for ii = 1:n
                w1 = kh(xb0 - xb(ii,:), hg); % (m,1)
                tmp = w1' .* gest(ii,:); % (1,m)
                tmp1(ii) = sum(tmp);
                w2 = kh(xb0 - xb(ii,:), hdg);
                tmp2(ii) = sum(w2);
            end
            gest_smooth(i,s) = sum(tmp1) / sum(tmp2);
        end
    end
    gest = gest_smooth;
end


end
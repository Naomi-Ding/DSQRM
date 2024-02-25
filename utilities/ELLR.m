%% ELLR (Empirical Log Likelihood Ratio) for beta
% function [Ui(beta)], i=1:n

function [ELLR, Uest] = ELLR(x, ally, beta, m, t, s, hx, hy, h, tau, G0, g_method)
% beta: (p,1)
% x: (n,p)
% ally: (n,m)

[n,p] = size(x);
xb = x*beta; % x.T * beta(sm), (n,1)

%% (1) Estimate g(X.T*Beta(s)) and g_dot(X.T*Beta(s))
% Eq(2.9)
y = ally(:,s); % (n,1), Y(s) at grid s
Gest = zeros(n,2);
if isequal(g_method, 'fminsearch')
    for i = 1:n
        xb0 = xb(i);
        % G0 = [g0, dg0]; % starting point for [g, dg]
        funG = @(G) polyG(G, xb, xb0, y, hy, tau);
        options = optimset('MaxFunEvals',10000, 'MaxIter',10000, 'Display', 'none'); %,'PlotFcns', @optimplotfval); % 'Display','Iter');
        Gest(i,:) = fminsearch(funG, G0(i,:), options); % estimators g and dg for object i
    end
    
elseif isequal(g_method, 'estimator')
    Gest(:,1) = G0(:,1);
    Gest(:,2) = G0(:,2);    
end

gest = Gest(:,1); % (n,1), g hat
dgest = Gest(:,2); % (n,1), gradient of g hat


mx = zeros(n,p);
U = zeros(n,p);
for i = 1:n % for object i
    %% (2) Estimate E{Xi|Xi.T*Beta(s)}
    % Eq(2.8)
    % --- locallinear --- %
    xb0 = xb(i);
    mxi = locallinear0(p, xb0, hx, xb, x); % E[xi|xi.T*beta(sm)], (1,p)
    mx(i, :) = x(i, :) - mxi; % xi-E[xi|xi.T*beta(sm)], (1,p)
    
    %% (3) Plug into Seff(Beta(s);Yi(sm),Xi} and Ui{Beta(s)}
    % Eq(2.5)
    delta = zeros(m,1);
    Seff = zeros(m,p); % Seff[beta(s);Yi(sm),Xi] as row vectors
    
    for ss = 1:m % at certain grid sm
        y = ally(:,ss); % (n,1) vector, Y(sm)
        delta(ss) = t(ss) - t(s); % scalar, (sm-s)
        Seff(ss,:) = mx(i,:) .* phi(y(i) - gest(i), tau) .* dgest(i).*kh(delta(ss),h); %(1,p) vector, Seff[beta(s);Yi(sm),Xi]
    end
    U(i,:) = mean(Seff,1); % mean of Seff[beta(s);Yi(sm),Xi] wrt m
    
end

ELLR = -2 * elmmatlab(U); % scalar, ELLR of beta(s)

if nargout > 1
    % disp('loading U')
    Uest = U;
end

end

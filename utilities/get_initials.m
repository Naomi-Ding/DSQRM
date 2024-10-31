function [beta0, g0, dg0] = get_initials(x, ally, t, hx, hy, h, h1, SIVC)

if ~exist('SIVC', 'var')
    SIVC = false;
end

[n,p] = size(x);
m = size(ally, 2);

%% initials for beta(s)
beta0 = zeros(p,m);
nd = 0;
for s = 1:m
    y = ally(:, s);
    tmp = Mimd1(x, y, nd);
    tmp = tmp / norm(tmp);
    beta0(:,s) = tmp;
end
beta0_smt = zeros(p,m);
for s = 1:m
    beta0_smt(:,s )= locallinear0(p, t(s), 0.05, t', beta0');
end
beta1_c = min(beta0_smt(1,:));
beta0_smt2 = beta0_smt + sign(beta1_c)*beta1_c + 1e-4;
beta0 = beta0_smt2 ./ sqrt(sum(beta0_smt2.^2, 1));
% beta0 = sign(beta0(1, 1)) * beta0;

if SIVC
    beta_sivc = zeros(p,m);
    for s = 1:m
        beta_sivc(:,s) = getBeta(x, ally, m, t, s, hx, hy, h, beta0(:,s));
    end
    beta_sivc_smt = zeros(p,m);
    for s = 1:m
        beta_sivc_smt(:,s) = locallinear0(p, t(s), 0.05, t', beta_sivc');
    end
    beta1_c = min(beta_sivc_smt(1,:));
    beta_sivc_smt2 = beta_sivc_smt + sign(beta1_c)*beta1_c + 1e-4;
    beta_sivc_smt2 = beta_sivc_smt2 ./ sqrt(sum(beta_sivc_smt2.^2, 1));
    beta0 = beta_sivc_smt2;
end

%% initials for g(xb) & dg(xb)
g0 = zeros(n,m);
dg0 = zeros(n,m);
allxb0 = x * beta0;
for i = 1:n
    for s = 1:m
        xb0 = x(i,:) * beta0(:,s);
        [g0(i,s), dg0(i,s)] = locallinear1(1, xb0, h1, allxb0(:), ally(:));
    end
end


end

%% ------- Global test statistic T ------- %%

function [Tn, feta_0, all_A, all_inv_meanB] = global_Tn(C, c, betaest, dgest, x, etaest, t, tau, hx, h3)

[n,m] = size(dgest);
[p,~] = size(betaest);

T = zeros(m,1);
d = sqrt(n) * (C * betaest - c); % (r,m)

%% (1) Estimate f_eta(0;s), the density of eta at x=0 for each s
feta_0 = eta_pdf(0, etaest, h3); % (1,m), density of eta at x=0 at each grid sm

allxb = x * betaest; %(n,m)
all_A = zeros(n,p,m);
all_inv_meanB = zeros(p,p,m);
for s = 1:m
    %% (2) Estimate Xi - E{Xi|Xi.T*Beta(s)}
    xb = allxb(:,s); %(n,1)
    mx = zeros(n,p);
    for i = 1:n
        xb0 = xb(i);
        mxi = locallinear0(p, xb0, hx, xb, x); % E[xi|xi.T*beta(sm)], (1,p)
        mx(i, :) = x(i, :) - mxi; % xi-E[xi|xi.T*beta(sm)], (1,p)
    end
    
    %% (3) Compute A, B for each s
    A = mx .* dgest(:,s); %(n,p)
    tmp1 = A' * A; % (p,p), sum(Ai*Ai')
    meanA = tmp1 ./ n; % (p,p), mean of Ai*Ai'
    tmp2 = mx' * (x .* (dgest(:,s).^2)); % (p,p), sum(Bi)
    meanB = tmp2 ./ n; % (p,p), mean of Bi
    inv_meanB = meanB^(-1);
    
    %% (4) Compute XI(s,s)
    XI = tau * (1-tau) * (feta_0(s)^(-2)) * inv_meanB * meanA * inv_meanB'; % (p,p), XI(s,s)
    tmp3 = C * XI * C'; % (r,r)
    T(s) = d(:,s)' * (tmp3^(-1)) * d(:,s); % scalar
    
    %% save A, meanB for each s
    all_A(:,:,s) = A;
    all_inv_meanB(:,:,s) = inv_meanB;
end

%% (5) Compute Tn
Tn = sum(T) * mean(diff(t)); % Global test statistic Tn


end
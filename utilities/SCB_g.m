%% SCB for index function g(z)

function Cg = SCB_g(x, ally, betaest, gest, dgest, tau, h1, zeta)

[n,m] = size(ally);
res = ally - gest;
yboot = gest + zeta .* res;

%% Compute [gest_boot, dgest_boot] with same betaest
gest_boot = zeros(n,m);
dgest_boot = zeros(n,m);
allxb = x * betaest; % (n,m)
xb = allxb(:); % (n*m,1) vector
options = optimset('MaxFunEvals',100, 'MaxIter',100, 'Display', 'none');
for i = 1:n % for object i
    for s = 1:m % at grid s
        xb0 = allxb(i,s);
        G0 = [gest(:,s), dgest(:,s)]; % (n,2), starting point of (g, dg)
        funG = @(G) minG(G, xb, xb0, yboot, m, h1, tau);
        G_est = fminsearch(funG, G0(i,:), options);
        gest_boot(i,s) = G_est(1);
        dgest_boot(i,s) = G_est(2);
    end
end

Gg = gest - gest_boot; % (n,m)
Cg = max(max(Gg)); % find Cg for each r

end
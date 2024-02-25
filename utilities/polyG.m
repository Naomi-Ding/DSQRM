%% Object function for estimators g and derivative of g %%

%% excluding sum for m
function funG = polyG(G, xb, xb0, y, hy, tau)
% G: 1 x 2
% g: scalar, g for object i
% dg: scalar, gradient of g for object i
% [g, dg] = G;
% xb = x*beta; % (n,1)
% y: (n,1), Y(s) at grid s

u = y - G(1) - G(2) * (xb - xb0); % (n,1)
w = kh(xb - xb0, hy); % (n,1)
f = rho(u, tau) .* w; % (n,1)
funG = sum(f);
end
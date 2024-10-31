%% Object Function for estimator g, eq(2.15)

function G = minG(G, xb, xb0, ally, m, h1, tau);

y = ally(:); % (n*m,1) vector

u = y - G(1) - G(2) * (xb - xb0); % (n*m,1) vector
w = kh(xb - xb0, h1); % (n*m,1) vector
funG = rho(u, tau) .* w; % (n,1) vector
G = sum(funG);

end

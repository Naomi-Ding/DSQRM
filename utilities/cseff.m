% Score function 

% INPUT
% x : n * p matrix
% y :  n * 1 vector
% b, hx, hy : scalar
% t : 1 * m vevtor
% beta : p * 1 vector

% OUTPUT
% yout : n * p matrix

% 02/25/2015
% Xinchao Luo

function [yout] = cseff(x, y, beta, hx, hy)

[n, p] = size(x); 
xb = x*beta; % x.T * beta(sm), (n,1)
mx0 = x;
yout = zeros(n, p);

gest = y;
for i=1:n % for object i
    xb0 = x(i, :)*beta; % xi.T * beta(sm), scalar
    gest(i) = locallinear0(1, xb0, hy, xb(:), y(:));  %g(xi.T*beta(sm))

end
vareps = y - gest; % epsilon_star(sm), (n,1)

for i=1:n
    x0 = xb(i);
    mx = locallinear0(p, x0, hx, xb, x); % E[xi|xi.T*beta(sm)], (1,p)
    mx0(i, :) = x(i, :) - mx; % xi-E[xi|xi.T*beta(sm)], (1,p)
end

fb = dlogeta(xb, y, hy); % dot_g(xi,T*beta(sm)), at grid sm
                         % derivative of local polynomial, (1,n)
% for i=1:n
%     x0 = xb(i);
%     tmp = locallinear0(1, x0, hx, xb, fb);
%     mfb(i) = fb(i) - tmp;
% end
for i=1:n
    % Seff[beta(sm);xi,Yi(sm)],eq(6), (1,p)
    yout(i,:) = vareps(i)*mx0(i,:)*fb(i); 
end

end
  

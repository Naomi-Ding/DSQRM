% Locallinear for beta0 

% INPUT
% p, x0, hx : scalar
% x : n-dimensional vector
% p, hx: scalar
% y : n * p matrix

% OUTPUT
% f : p-dimensional vector

% 02/25/2015
% Xinchao Luo

function [f] = locallinear0(p, x0, hx, x, y)

eps = 1e-8;
u = x - x0; % (n,1) vector
w = kh(u, hx); % (n,1) vector
tmp1 = w.*u; 
tmp2 = tmp1.*u;
a = sum(w);
b = sum(tmp1);
d = sum(tmp2);
deno = a*d-b^2;
f = zeros(1, p);
for j=1:p
    tmp1 = w.*y(:, j);
    tmp2 = tmp1.*u;
    t1 = sum(tmp1);
    t2 = sum(tmp2);
    num1 = d*t1-b*t2;
    f(j) = num1/(deno+eps);
end

end


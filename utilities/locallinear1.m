% Locallinear for beta0 and beta1 

% INPUT
% p, x0, hx : scalar
% x : n-dimensional vector
% p, hx: scalar
% y : n * p vector

% OUTPUT
% f1, f2 : p-dimensional vector

% 02/25/2015
% Xinchao Luo

function [f1, f2] = locallinear1(p, x0, hy, x, y)

eps = 1e-8;
u = x - x0;
f1 = zeros(1, p);
f2 = f1;
w = kh(u, hy);
tmp1 = w.*u;
tmp2 = tmp1.*u;
a = sum(w);
b = sum(tmp1);
d = sum(tmp2);
deno = a*d-b^2;
for j=1:p
    tmp1 = w.*y(:, j);
    tmp2 = tmp1.*u;
    t1 = sum(tmp1);
    t2 = sum(tmp2);
    num1 = d*t1-b*t2;
    f1(j) = num1/(deno+eps);
    num2 = -b*t1+a*t2;
    f2(j) = num2/(deno+eps);
end 

end


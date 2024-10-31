% derivative of local polynomial

% INPUT
% hy : scalar
% xb : n-dimensional vector
% y : n * p vector

% OUTPUT
% f : n-dimensional vector

% 08/05/2015
% Xinchao Luo

function [f] = dlogeta(xb, y, hy)

n = length(y);
f = zeros(1, n);
for i=1:n
   x0 = xb(i);
   [~, f2] = locallinear1(1, x0, hy, xb, y);
   f(i) = f2; 
end

end
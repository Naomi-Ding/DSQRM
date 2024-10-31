% Sum of weighted sum of score function 

% INPUT
% x : n * p matrix
% ally :  n * m matrix
% s, m, b, hx, hy, h : scalar
% t : 1 * m vevtor
% beta : p * 1 vector

% OUTPUT
% fval : 1 * p vector

% 07/29/2015
% Xinchao Luo

function [fval] = sumseff(x, ally, m, t, s, beta, hx, hy, h)

yout = allcseff(x, ally, m, t, s, beta, hx, hy, h); % sum of M grids weighted Seff, (n,p)
fval = mean(yout, 1); % mean of each column, (1,p)
% find beta_hat(s)by solving sumSeff=0, equivalent to solving meanSeff=0
end
  

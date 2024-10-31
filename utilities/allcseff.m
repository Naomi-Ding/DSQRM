% Weighted sum of score function 

% INPUT
% x : n * p matrix
% ally :  n * m matrix
% s, m, hx, hy : scalar
% t : 1 * m vevtor
% beta : p * 1 vector

% OUTPUT
% yout : n * p vector

% 02/25/2015
% Xinchao Luo

function [yout] = allcseff(x, ally, m, t, s, beta, hx, hy, h)

delta = zeros(m, 1);
for ss=1:m
    delta(ss) = t(ss) - t(s);
end
kernel = kh(delta, h); %weight function, w(sm,s)=Kh(sm-s), (m,1)

yout = 0;
for ss=1:m % at some certain grid sm
    if (kernel(s)>0)
        y = ally(:, ss); % y(sm)
        tmp = cseff(x, y, beta, hx, hy); %Seff[beta(sm);X,Y(sm)],(n,p)
        yout = yout + kernel(ss)*tmp; % sum of weighted Seff at grid s1,...,sm
    end
end

end
  

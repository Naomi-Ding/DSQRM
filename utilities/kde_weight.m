%% weight function w(t,hf) for kernel density estimator
function w = kde_weight(t, kernel, hf)
% t: 1 x N points in [0,1]
% hf: bandwidth, hf < 1/2

% w: 1 x N, weight function for at each point ti

if isequal(kernel, 'ep')
%     K = @(x,vi) kh((x - vi)/hf, 1);
    K = @(u) 0.75 * (1-u.^2) .* (abs(u)<1);
elseif isequal(kernel, 'gaussian')
%     exp(-tmp(i)^2/2) / (2*pi * h);
    K = @(u) exp(-u.^2/2) / sqrt(2*pi);
end

N = length(t);
w = ones(1,N);
for i = 1:sum(t < hf)
    w(i) =  1 / (integral(K, -t(i)/hf, 1));
end
for i = (N - sum(t > t(end)-hf) + 1) : N
    w(i) = 1 / (integral(K, -1, (t(end)-t(i))/hf));
end

end

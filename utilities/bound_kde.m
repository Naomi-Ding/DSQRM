%% boundary-corrected kernel density estimator
function [f,t] = bound_kde(v, kernel, hf)
% v: n x N, voxel-wise imaging measurement {v_i(t_N)}
% kernel: choice of kernel function, 'gaussian' or 'ep'
% hf: scalar or n x 1 vector, bandwidth for each subject, hf<1/2

% f: n x N, density functions {f_i(t)}
% t: n x N, support for density, N points in [0, upper bound]

[n,N] = size(v);
% m = length(t);

if isequal(kernel, 'ep')
    %     K = @(x,vi) kh((x - vi)/hf, 1);
    K = @(u) 0.75 * (1-u.^2) .* (abs(u)<1);
elseif isequal(kernel, 'gaussian')
    %     exp(-tmp(i)^2/2) / (2*pi * h);
    K = @(u) exp(-u.^2/2) / sqrt(2*pi);
end

if length(hf) == 1 % universe bandwidth
    hf = repmat(hf, n, 1);
end

w = @(t, hf) kde_weight(t, kernel, hf); % 1 x N, weight function

bound_u = max(v,[],2) + 0.001; % n x 1, upper bound
bound_l = min(v,[],2) - 0.001; % n x 1, lower bound
f = zeros(n,N);
t = zeros(n,N);
for i = 1:n
    %     i
    %     ti = linspace(0, bound(i), N); % 1 x N, coordinates for density fi
    ti = linspace(bound_l(i), bound_u(i), N); % 1 x N, coordinates for density fi
    Kw = K((ti - v(i,:)')/hf(i)) .* w(ti, hf(i)); % N x N
    kw0 = trapz(ti, Kw, 2); % N x 1, integral w.r.t. t
    f(i,:) = sum(Kw, 1) / sum(kw0);
    t(i,:) = ti;
end

end

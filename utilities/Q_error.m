%% error-contaminated random quantile functions
function Q = Q_error(U, ratio, tau, x, lambda, eps)
% U: n x N in [0,1], sampled N points for n subjects from Uniform(0,1)
% tau: scalar in [0,1], tau_th quantile of eta(s)
% x: n x p, covariates
% lam: 1 x 2, lambda, sigma^2 of epsilon for generating eta(s)
% eps: n x 2, epsilon sampling from standard normal distribution 

% Q: n x N, quantile function Q(u|x) 

beta0 = @(s) [1+s.^2; ratio * (1-s).^2]; % p x N
beta = @(s) beta0(s) ./ vecnorm(beta0(s)); % p x N
g = @(u) sin(2*pi*u).*(u+0.5);

psi1 = @(s) cos(2*pi.*s');
psi2 = @(s) sin(pi.*s');

[n, N] = size(U);
sort_U = sort(U, 2); % sort elements of each subject (each row)

t = linspace(0,1,N); % 1 x N
eta_t = [lambda(1)*eps(:,1),lambda(2)*eps(:,2)]*[psi1(t) psi2(t)]'; % n x N
eta_t = eta_t - quantile(eta_t, tau,1); % n x N
g_t = g(x*beta(t)); % n x N
integ = exp(g_t + eta_t); % n x N

Q = zeros(n,N);
for i = 1:n 
	u = sort_U(i,:); % 1 x N
	integ_u = interp1(t, integ(i,:), u); % 1x N
	for k = 1:N 
		Q(i, k) = trapz([0, u(1:k)], [integ(i,1), integ_u(1:k)]); 
	end
end 

end 
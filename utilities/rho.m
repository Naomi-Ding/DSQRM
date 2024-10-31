%% Check Function %%

function [f] = rho(u, tau)

% f = sum(abs(u.*(tau - (u<0))));
f = u.* (tau - (u < 0));

end
%% Density function of eta(s)

function [f] = eta_pdf(x,etaest,h3)

[n,m] = size(etaest);
k = length(x); % length of vector x
f = zeros(k,m);

for s = 1:m
    kernel = zeros(k,n);
    for i = 1:n
        kernel(:,i) = kh((x - etaest(i,s)),h3); %(k,1), same shape with x
    end
    f(:,s) = mean(kernel,2); %(k,1)
end

end

    
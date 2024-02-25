%% SCV function for h3

function f = SCV(etaest, h3)

[n,m] = size(etaest);
I = zeros(m,1);

for s = 1:m
    tmp1 = zeros(n,1);
    for i = 1:n
        w = etaest(i,s) - etaest(:,s); % (n,1), each row is eta_i(s) - eta_j(s)
        kernel = kh(w, h3); % (n,1), each row is Kh3[eta_i(s) - eta_j(s)]
        kernel(i) = []; % (n-1,1), j does not equal to i
        tmp1(i) = sum(kernel, 1); % scalar, sum_j[Kh3(eta_i - eta_j)]
    end
    tmp2 = sum(tmp1,1); % scalar, sum_i
    I2 = tmp2 * 2 / (n * (n-1)); % estimation of I2(sm)
    
    f_eta = @(x) eta_pdf(x, etaest(:,s), h3); % f_eta(x;sm)
    f_eta2 = @(x) f_eta(x).^2; % [f_eta(x;sm)]^2
    I1 = integral(f_eta2, -Inf, Inf, 'ArrayValued', true);
    I(s) = I1 - I2; % (m,1), I1-I2 at grid sm
end

f = sum(I); % object function

end

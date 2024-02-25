function [betaest, gest] = smooth_estimators(betaest, gest, t, x, hb, hg, hdg)
if ~exist('hb')
    hb = 0.05;
end
if ~exist('hg')
    hg = 0.05;
end
if ~exist('hdg')
    hdg = 0.05;
end

[p,m] = size(betaest);
n = size(x, 1);

betaest_smooth = zeros(p,m);
for s = 1:m
    betaest_smooth(:,s) = locallinear0(p, t(s), hb, t', betaest');
end
betaest_smooth = betaest_smooth ./ sqrt(sum(betaest_smooth.^2,1));
betaest = betaest_smooth;

gest_smooth = zeros(n,m);
xb = x * betaest;
for i = 1:n
    for s = 1:m
        xb0 = xb(i,s); % scalar
        tmp1 = zeros(n,1);
        tmp2 = zeros(n,1);
        for ii = 1:n
            w1 = kh(xb0 - xb(ii,:), hg); % (m,1)
            tmp = w1' .* gest(ii,:); % (1,m)
            tmp1(ii) = sum(tmp);
            w2 = kh(xb0 - xb(ii,:), hdg);
            tmp2(ii) = sum(w2);
        end
        gest_smooth(i,s) = sum(tmp1) / sum(tmp2);
    end
end
gest = gest_smooth;

end
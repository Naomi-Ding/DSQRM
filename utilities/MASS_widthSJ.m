%% Functions modified from R library MASS::width.SJ
function [h] = MASS_widthSJ(x, nb, method)
if nargin == 1
    nb = 1000;
    method = "ste";
end

SDh = @(x, h, n, d) VR_phi4_bin(n, length(x), d, x, h);
TDh = @(x, h, n, d) VR_phi6_bin(n, length(x), d, x, h);

n = length(x);
[cnt, d] = VR_den_bin(n, nb, x);
hmax = 1.144 * sqrt(var(x)) * n^(-1/5);
lower = 0.1 * hmax;
upper = hmax;
scale = min(sqrt(var(x)), iqr(x) / 1.349);
a = 1.24 * scale * n^(-1/7);
b = 1.23 * scale * n^(-1/9);
c1 = 1 / (2 * sqrt(pi) * n);
TD = TDh(cnt, b, n, d);
alph2 = 1.357 * (SDh(cnt, a, n, d) / TD)^(1/7);

if isequal(method, "dpi")
    res = (c1 / SDh(cnt, (2.394 / (n * TD))^(1/7), n, d))^(1/5);
else
    fSD = @(h, x, alph2, c1, n, d) (c1 / SDh(x, alph2 * h^(5/7), n, d))^(1/5) - h;
    itry = 1;
    if fSD(lower, cnt, alph2, c1, n, d) * fSD(upper, cnt, alph2, c1, n, d) > 0
        if itry > 99
            error("no solution in the specified range of bandwidths");
        end
        if mod(itry, 2)
            upper = upper * 1.2;
        else
            lower = lower / 1.2;
        end
        itry = itry + 1;
        % error("no solution in the specified range of bandwidths");
    end
    options = optimset('TolX',0.1*lower); %, "Display", "iter");
    res = fminbnd(@(h) abs(fSD(h, cnt, alph2, c1, n, d)), lower, upper, options);
end

h = 4 * res;
end


function [u] = VR_phi4_bin(n, nb, d, x, h)
DELMAX = 1000;
nn = n;
nbin = nb;
sum = 0;
for i = 1:nbin
    delta = i * d / h;
    delta = delta * delta;
    if delta >= DELMAX
        break;
    end
    term = exp(-delta/2) * (delta * delta - 6 * delta + 3);
    sum = sum + term * x(i);
end
sum = 2 * sum + nn * 3;
u = sum / (nn * (nn - 1) * (h^5) * sqrt(2 * pi));
end


function [u] = VR_phi6_bin(n, nb, d, x, h)
DELMAX = 1000;
nn = n;
nbin = nb;
sum = 0;
for i = 1:nbin
    delta = i * d / h ;
    delta = delta * delta;
    if delta > DELMAX
        break;
    end
    term = exp(- delta / 2) * (delta^3 - 15 * (delta^2) + 45 * delta -15);
    sum = sum + term * x(i);
end
sum = 2 * sum - 15*nn;
u = sum / (nn * (nn-1) * (h^7) * sqrt(2 * pi));
end




function [cnt, dd] = VR_den_bin(n, nb, x)
INT_MAX = 1e9;
nn = n;
cnt = zeros(1, nb);
xmin = min(x);
xmax = max(x);
% for i = 1:nb
%     cnt(i) = 0;
% end
% xmin = xmax = x(1);
% for i = 2:nn
%     xmin = min(xmin, x(i));
%     xmax = max(xmax, x(i));
% end
rang = (xmax - xmin) * 1.01;
dd = rang / nb;
for i = 2:nn
    ii = floor(x(i) / dd);
    for j = 1:(i-1)
        jj = floor(x(j) / dd);
        iij = abs(ii - jj) + 1;
        if cnt(iij)== INT_MAX
            error("maximum count exceeded in pairwise distance binning");
        end
        cnt(iij) = cnt(iij) + 1;
    end
end
end
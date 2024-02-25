% CV for h1 

% INPUT
% x : n * p matrix
% betaest : p * m matrix
% allxb, ally : n * m matrix

% OUTPUT
% h1 : scalar

% 08/25/2015
% Xinchao Luo

function [h1] = cvh1(x, betaest, allxb, ally)

[n, m] = size(allxb);
% srange = range(allxb(:));
srange = max(allxb(:) - min(allxb(:)));
nh = 20; CV = zeros(nh,1);
hmin = srange/m; 
hmax = srange/20; 
vh = logspace(log10(hmin),log10(hmax),nh);
gest = zeros(n, m);

for ii=1:nh
    h = vh(ii);
    ncv = 0;
    for i=1:n
        for s=1:m
            ncv = ncv + 1;
            xb0 = x(i, :)*betaest(:, s);
            xb_cv = allxb(:); xb_cv(ncv) = [];
            y_cv = ally(:); y_cv(ncv) = [];
            gest(i, s) = locallinear0(1, xb0, h, xb_cv, y_cv);        
        end
    end
    CV(ii) = mean(mean((ally - gest).^2));
end

[~,flag] = min(CV);
h1 = vh(flag);

end


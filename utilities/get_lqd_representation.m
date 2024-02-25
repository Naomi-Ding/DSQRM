function [lqd, lqdSup, hf, fhat, supp] = get_lqd_representation(v, m)
if nargin == 1
    m = 100;
end

if ~iscell(v)
    n = size(v, 1);
    %% Density Estimation & Obtain LQD
    %% (1) Bandwidth Selection:
    hf = zeros(1, n);
    for i = 1:n
        hf(i) = MASS_widthSJ(v(i, :));
    end
    
    %% (2) Boundary-corrected Kernel Density Estimator
    [f,supp] = bound_kde(v, 'ep', hf);
    fhat = f;
    fhat(f==0) =  f(f==0) + 0.0001;
    
    %% (3) Obtain LQD for fhat
    lqdSup = linspace(0, 1, m);
    lqd = zeros(n,m);
    for i = 1:n
        lqd(i, :) = dens2lqd(fhat(i,:), supp(i,:), m);
    end
    
else
    n = length(v);
    v_cell = v;
    %% Density Estimation & Obtain LQD
    %% (1) Bandwidth Selection & (2) Boundary-corrected Kernel Density Estimator
    hf = zeros(1, n);
    f0 = cell(n,1);
    supp0 = cell(n,1);
    maximum = 0;
    for i = 1:n
        vi = double(v_cell{i});
        %% (1) Bandwidth Selection
        hf(i) = MASS_widthSJ(vi);
        %% (2) Boundary-corrected Kernel Density Estimator
        [fi, suppi] = bound_kde(vi', 'ep', hf(i));
        fi(fi==0) = fi(fi==0) + 0.0001;
        suppi = suppi - min(suppi);
        maximum = max(max(suppi), maximum);
        f0{i} = fi;
        supp0{i} = suppi;
    end
    fhat = cell(n,1);
    supp = cell(n,1);
    lqdSup = linspace(0, 1, m);
    lqd = zeros(n,m);
    for i = 1:n
        suppi = supp0{i} / maximum; % normalize the support
        fi = f0{i};
        c = trapz(suppi, fi);
        fhat{i} = fi / c;
        supp{i} = suppi;
        lqd(i, :) = dens2lqd(fhat{i}, suppi, m);
    end
end

end

%% Simulation results
function [T_all, all_betaest, all_gest, all_dgest, all_gest_inv, ...
    all_betaest_SIVC, all_gest_SIVC, all_gest_inv_SIVC] = ...
    simu_main(n, p, m, N, nsimu, tau_set, SIVC)
addpath('./utilities')
if ~exist('simu_results', "dir")
    mkdir('simu_results');
end

tic;
%% Data generation: sampling from quantile functions with AR(1) correlation
disp("Generate Simulation Dataset");
[simu_x, simu_g, all_v, beta, g, dg, t, hx, hy, h] = simu_data_gen(n, p, m, N, nsimu);
beta0 = beta(t);
toc;
% save the dataset
dataname = sprintf("simu_results/simu_n%d_N%d_m%d_p%d_nsimu%d.mat", n, N, m, p, nsimu);
save(dataname, "simu_x", "simu_g", "all_v", "beta0", "t", "hx", "hy", "h");

%% Density Estimation & LQD Transformation & Model fitting
disp("Extract density estimators & LQD representations & Estimate parameters");
all_fhat = zeros(n,N, nsimu);
all_fsupp = zeros(n,N, nsimu);
all_hf = zeros(nsimu,n);
simu_ally = zeros(n,m, nsimu);
ntau = length(tau_set);
all_betaest = zeros(p,m,nsimu, ntau);
all_gest = zeros(n,m,nsimu, ntau);
all_dgest = zeros(n,m,nsimu, ntau);
all_gest_inv = zeros(n, m,nsimu, ntau);
parfor nn = 1:nsimu
    disp(compose("nn = %d\n", nn));
    x = simu_x(:,:,nn);
    v = all_v(:,:,nn);
    g0 = g(x*beta0);
    dg0 = dg(x*beta0);
    
    [fhat, f_support, hf, ally, betaest, gest, dgest, gest_inv] = ...
        DSQRM(x, v, m, tau_set, hx, hy, h, [], beta0, g0, dg0, false, 0);
    all_fhat(:,:,nn) = fhat;
    all_fsupp(:,:,nn) = f_support;
    all_hf(nn,:) = hf;
    simu_ally(:,:,nn) = ally;
    
    all_betaest(:,:,nn,:) = betaest;
    all_gest(:,:,nn,:) = gest;
    all_dgest(:,:,nn,:) = dgest;
    all_gest_inv(:,:,nn,:) = gest_inv;
end
toc;
save(dataname, "all_betaest", "all_gest", "all_dgest", "all_gest_inv", ...
    "simu_ally", "all_fhat", "all_fsupp", "all_hf", "-append");

%% TRUE density f & the inverse of true link function
truef = zeros(n, N, nsimu);
g_inv = zeros(n,m,nsimu);
parfor nn = 1:nsimu
    for i = 1:n
        truef(i, :, nn) = lqd2dens(simu_g(i, :,nn), t, all_fsupp(i,:,nn));
    end
    g_inv(:,:,nn) = MakeDENsample(simu_g(:,:,nn));
end
save(dataname, "truef","g_inv", "-append");

%% Calculate mean & std of ISE
tmp = (all_fhat - truef).^2; tmp = mean(tmp, 2);  tmp = mean(tmp, 1);%tmp = tmp.^0.5;
ISE_f = [mean(tmp), std(tmp)]';
tmp = (all_betaest - repmat(beta0, 1, 1, nsimu, ntau)).^2; tmp = squeeze(mean(tmp, 2));
ISE_beta = [mean(tmp,2), std(squeeze(tmp), 0,2)]; % (2,2,ntau)
tmp = (all_gest- repmat(simu_g, 1, 1, 1, ntau)).^2; tmp = mean(tmp, 2); tmp = squeeze(mean(tmp, 1));
ISE_g = [mean(tmp); std(tmp, 0, 1)]'; % (ntau, 2)
tmp = (all_gest_inv - repmat(g_inv, 1,1,1,ntau)).^2; tmp = mean(tmp, 2); tmp = squeeze(mean(tmp, 1));
ISE_g_inv = [mean(tmp); std(tmp)]'; % (ntau, 2)

varnames = {'n','N','ISE_f', 'm', 'tau', 'ISE_beta1', 'ISE_beta2', 'ISE_g', 'ISE_Psi^{-1}(g)'};
T = table(n*ones(ntau,1), N*ones(ntau,1), repmat(ISE_f', ntau,1), ...
    m*ones(ntau,1), tau_set', reshape(squeeze(ISE_beta(1,:,:)), 2, [])', ...
    reshape(squeeze(ISE_beta(2,:,:)), 2, [])', ISE_g, ISE_g_inv,  'VariableNames', varnames);
writetable(T, "simu_results/simu_estimation_error.csv", "WriteMode", "append");



%% Comparison: SIVC
if SIVC
    disp("Estimation by SIVC");
    all_betaest_SIVC = zeros(p,m,nsimu);
    all_gest_SIVC = zeros(n,m,nsimu);
    all_gest_inv_SIVC = zeros(n,m,nsimu);
    parfor nn = 1:nsimu
        disp(compose("nn = %d\n", nn));
        x = simu_x(:,:,nn);
        ally = simu_ally(:,:,nn);
        [betaest, gest] = get_SIVC_estimate(x, ally, t, hx, hy, h, beta0);
        % Get the inverse of estimated link
        all_gest_inv_SIVC(:,:,nn) = MakeDENsample(gest);
        all_betaest_SIVC(:,:, nn) = betaest;
        all_gest_SIVC(:,:,nn) = gest;
    end
    toc;
    
    % mean & std of ISE
    tmp = (all_betaest_SIVC - repmat(beta0, 1, 1, nsimu)).^2; tmp = mean(tmp, 2);
    ISE_beta_SIVC = [mean(tmp,3), std(squeeze(tmp), 0,2)]; % (2,2)
    tmp = (all_gest_SIVC- simu_g).^2; tmp = mean(tmp, 2); tmp = mean(tmp, 1);
    ISE_g_SIVC = [mean(tmp), std(squeeze(tmp), 0, 1)]; % (1, 2)
    tmp = (all_gest_inv_SIVC - g_inv).^2; tmp = mean(tmp, 2); tmp = mean(tmp, 1);
    ISE_g_inv_SIVC = [mean(tmp), std(squeeze(tmp), 0, 1)]; % (1, 2)
    
    % save the results & display the errors
    save(dataname, 'all_betaest_SIVC', 'all_gest_SIVC', 'all_gest_inv_SIVC', '-append');
    varnames = {'n','N','ISE_f', 'm', 'tau', 'ISE_beta1', 'ISE_beta2', 'ISE_g', 'ISE_Psi^{-1}(g)'};
    T_SIVC = table(n, N, ISE_f', m, "SIVC", ISE_beta_SIVC(1,:), ISE_beta_SIVC(2,:), ...
        ISE_g_SIVC, ISE_g_inv_SIVC,  'VariableNames', varnames);
    writetable(T_SIVC, "simu_results/simu_estimation_error.csv", "WriteMode", "append");
    T_all = [T; T_SIVC];
    disp(T_all);
else
    T_all = T;
    all_betaest_SIVC = [];
    all_gest_SIVC = [];
    all_gest_inv_SIVC = [];
end

% %% save the results
% dataname = sprintf("simu_results/simu_n%d_N%d_m%d_p%d_nsimu%d.mat", n, N, m, p, nsimu);
% save(dataname, "simu_x", "simu_ally", "hx", "hy", "h", "t", "simu_g", ...
%     "truef","g_inv", "beta0", "ISE_f");
end

%% table 2: partial results
clear all;
addpath('./utilities')

%%% settings
n = 100;
p = 2;
m = 100;
N = 500;
nsimu = 200; % only generate 200 datasets to replicate the n-th set
% tau_set = 0.1:0.2:0.9;
tau_set = 0.5;
SIVC = false; 

[T_all, all_betaest, all_gest, all_dgest, all_gest_inv, ...
	all_betaest_SIVC, all_gest_SIVC, all_gest_inv_SIVC] = ...
	simu_main(n, p, m, N, nsimu, tau_set, SIVC);

disp(T_all);
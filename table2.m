%% table 2: Estimation results from 200 simulated datasets for N = 500
clear all;
addpath('./utilities')
N = 500;
p = 2;
nsimu = 200;
SIVC = true; % if true: run SIVC on the same dataset as well

%% settings: n=100, m=100
n = 100; % sample size
m = 100; % number of grids for lqd functions
tau_set = 0.5;
% tau_set = 0.1:0.2:0.9; % all
[T_all, all_betaest, all_gest, all_dgest, all_gest_inv, ...
	all_betaest_SIVC, all_gest_SIVC, all_gest_inv_SIVC] = ...
	simu_main(n, p, m, N, nsimu, tau_set, SIVC);
disp(T_all);
pause;

%% settings: n=100, m=200
n = 100;
m = 200;
tau_set = 0.5;
% tau_set = 0.1:0.2:0.9; % all
[T_all, all_betaest, all_gest, all_dgest, all_gest_inv, ...
	all_betaest_SIVC, all_gest_SIVC, all_gest_inv_SIVC] = ...
	simu_main(n, p, m, N, nsimu, tau_set, SIVC);
disp(T_all);
pause;

%% settings: n=200, m=100
n = 200;
m = 100;
tau_set = 0.5;
% tau_set = 0.1:0.2:0.9; % all
[T_all, all_betaest, all_gest, all_dgest, all_gest_inv, ...
	all_betaest_SIVC, all_gest_SIVC, all_gest_inv_SIVC] = ...
	simu_main(n, p, m, N, nsimu, tau_set, SIVC);
disp(T_all);
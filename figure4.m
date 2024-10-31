%% GBM dataset, tau = 0.5, estimation & SCB
clear all;
%% settings
m = 100;
tau = 0.5;
smooth = true;
R = 200;
a = 0.05;
data_folder = "Software_DSQRM/examples/";
addpath("utilities/");

%% input: images & covariates
tic;
disp("Step 1. Extract density estimators & LQD representations");
% load the covariates
df = readtable(fullfile(data_folder, "TCGA_GBM_covariates.csv"));
sample_name = df.(1);
x0 = table2array(df(:, 2:end));
x0name = df.Properties.VariableNames(2:end);

% load the image
image_folder = fullfile(data_folder, "TCGA_flair_single_slice/");
[tumor_pixels, all_tumor_coord, sub_tumor_ratios, all_images, ...
    sub_tumor_coord, rm_idx] = extract_pixels_from_NIFTI(sample_name, image_folder);
% extract density & LQD
[lqd, lqdSup, hf, fhat, supp] = get_lqd_representation(tumor_pixels, m);
toc;

% prepare the design matrix
% normalize covariates to be in the range [0,1]
xidx = max(x0, [], 1) > 1;
x0(:, xidx) = normalize(x0(:, xidx), 'range');
if ~isempty(rm_idx)
    x0(rm_idx, :) = [];
end
p0 = size(x0, 2);
ptumor = size(sub_tumor_ratios, 2);
xdesign = [sub_tumor_ratios(:, end-1), x0, sub_tumor_ratios(:, 1:(end-2))];
% xname = [["subtype", num2str(ptumor - 1)], x0name];
xname = ["edema", x0name, "necrosis"];
[n, p] = size(xdesign);

%% bandwidth & initials
std_x = std(xdesign(:));
hx = std_x * n^(-1/3);
hy = std_x * n^(-1/5);
h_beta = (m)^(-1/5) * 0.7;
h3 = 0.4355;
h_g = 0.15;
h_eta = 0.02;

init = load(fullfile(data_folder, "initials.mat"));
beta0 = init.beta0;
g0 = init.g0;
dg0 = init.dg0;

%% estimate & SCB
disp("Step 2. Estimate the functional coefficients & link function");
[betaest, gest, dgest, ~] = get_estimate(xdesign, lqd, lqdSup, tau, hx, hy, ...
    h_beta, beta0, g0, dg0, smooth, 1, h_g);

disp("Step 3. Confidence Bands for the functional coefficients & link function");
ystar = lqd - gest;
[SCB_beta] = inference_beta(xdesign, ystar, betaest, dgest, lqdSup, tau, hx, R, a, h_eta, 100*h3);
SCB_g = 0.2001;
inf_data = load(fullfile(data_folder, "results.mat"));
p_vals = inf_data.p_vals;
toc;

gest_smooth = zeros(n,m);
hg = 0.05;
hdg = hg;
xb = xdesign * betaest;
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


%% plot the results
betaL = betaest - SCB_beta;
betaU = betaest + SCB_beta;
gL = gest - SCB_g;
gU = gest + SCB_g;

figure(4); clf;
sgtitle("$\tau=0.5$", 'Interpreter', 'latex', 'FontSize',18);
for k = 1:p
    ax_k = subplot(2, ceil(p/2), k);
    plot(ax_k, lqdSup, betaest(k,:), 'k-', 'LineWidth', 2);
    hold(ax_k, 'on');
    plot(ax_k, lqdSup, betaL(k,:), 'k--', 'LineWidth', 2);
    plot(ax_k, lqdSup, betaU(k,:), 'k--', 'LineWidth', 2);
    txt = [xname{k}, '(', num2str(p_vals(k)), ')'];
    if p_vals(k) < a
        txt_color = 'r';
    else
        txt_color = 'k';
    end
    title(ax_k, txt, "FontSize", 14, 'Color', txt_color, 'Interpreter', 'none');
    hold(ax_k, 'off');
end
ax_g = subplot(2, ceil(p/2), p+1);
allxb = xdesign * betaest;
xb_plot = allxb(:); [xb_plot, index] = sort(xb_plot);
g_plot = gest(:); g_plot = g_plot(index);
gL_plot = gL(:); gL_plot = gL_plot(index);
gU_plot = gU(:); gU_plot = gU_plot(index);
plot(ax_g, xb_plot, g_plot, "k-", "LineWidth", 2);
hold on;
plot(ax_g, xb_plot, gL_plot, 'k--', xb_plot, gU_plot, 'k--', 'LineWidth',2);
xlabel(ax_g, 'x^T\beta(s)', 'FontSize', 14); title(ax_g, 'g(x^T\beta(s))', 'FontSize', 14);
hold(ax_g, 'off');
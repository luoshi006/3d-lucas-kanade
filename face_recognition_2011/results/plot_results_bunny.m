% plot_results_AR
% reproduces Fig 5a of the paper 
% G. Tzimiropoulos, S. Zafeiriou and M. Pantic, "Robust and Efficient Parametric Face Alignment", ICCV 2011

clear; clc;
n_freq_tests = 1;
load results10.mat; %load results
results = mean(results, 2);
results = squeeze(results);

%%
affine_ECC_ic = results(:, 1) / n_freq_tests;
affine_GradientCorr_Euclidean_Split_ic_3d = results(:, 2) / n_freq_tests;
affine_GradientCorr_Euclidean_ic_3d = results(:, 3) / n_freq_tests;
affine_GradientImages_Split_ic_3d = results(:, 4) / n_freq_tests;
affine_GradientImages_ic_3d = results(:, 5) / n_freq_tests;
affine_ic = results(:, 6) / n_freq_tests;
affine_ic_irls = results(:, 7) / n_freq_tests;

var = 1:10;
figure; 
plot(var, affine_ECC_ic, 'black--diamond', ...
     var, affine_GradientCorr_Euclidean_Split_ic_3d, 'yellow:^',  ...
     var, affine_GradientCorr_Euclidean_ic_3d, 'c-.x', ...
     var, affine_GradientImages_Split_ic_3d, 'r:*', ...
     var, affine_GradientImages_ic_3d, 'green:^', ...
     var, affine_ic, 'blue-s', ...
     var, affine_ic_irls, 'red-s', ...
     'MarkerSize',11, 'linewidth', 2); 
grid on
set(gca, 'FontSize', 0.1)
set(gca, 'FontWeight', 'bold')
xtick = var; 
ytick = 0:0.1:1;
set(gca, 'xtick', xtick);
set(gca, 'ytick', ytick);
ylabel('Frequency of Convergence', 'Interpreter','tex', 'fontsize', 15);
xlabel('Point Standard Deviation', 'Interpreter','tex', 'fontsize', 15);
legend(gca, 'ECC-ic', 'GradCorrSplit-ic',  'GradCorr-ic', 'GradImSplit-ic', 'GradIm-ic', 'LK-ic', 'IRLS-ic');
title('Brain with Stanford Bunny', 'Interpreter', 'tex', 'fontsize', 15);

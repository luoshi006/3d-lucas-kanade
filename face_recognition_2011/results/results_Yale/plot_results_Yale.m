% plot_results_AR
% reproduces Fig 5a of the paper 
% G. Tzimiropoulos, S. Zafeiriou and M. Pantic, "Robust and Efficient Parametric Face Alignment", ICCV 2011

clear; clc; close all

load resultsYaleSmoothing.mat; %load results

results = mean(results,2);
results = mean(results,1);
results = squeeze(results);

affine_ic = results(:,1)/100;
affine_ECC_ic = results(:,2)/100;
affine_ic_irls = results(:,3)/100;
affine_GaborFourier_ic = results(:,4)/100;
affine_GradientImages_ic = results(:,5)/100;
affine_GradientCorr_ic = results(:,6)/100;

var = 1:10;
figure; plot(var, affine_ic, 'black--diamond', var, affine_ECC_ic, 'yellow:^',  var, affine_ic_irls, 'c-.x', var, affine_GaborFourier_ic, 'r:*', var, affine_GradientImages_ic, 'green:^', var, affine_GradientCorr_ic, 'blue-s',  'MarkerSize',11, 'linewidth', 2); grid on
set(gca, 'FontSize', 0.1)
set(gca, 'FontWeight', 'bold')
xtick = var; 
ytick = 0:0.1:1;
set(gca, 'xtick', xtick);
set(gca,'ytick', ytick);
ylabel('Frequency of Convergence', 'Interpreter','tex', 'fontsize', 15)
xlabel('Point Standard Deviation', 'Interpreter','tex', 'fontsize', 15)
legend('LK-ic', 'ECC-ic',  'IRLS-ic', 'GaborFourier-ic', 'GradientImages-ic', 'GradientCorr-ic')
title('Yale: with Smoothing','Interpreter','tex', 'fontsize', 15)

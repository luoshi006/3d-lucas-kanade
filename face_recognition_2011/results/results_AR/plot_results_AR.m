% plot_results_AR
% reproduces Fig 5b and Fig 5c of the paper 
% G. Tzimiropoulos, S. Zafeiriou and M. Pantic, "Robust and Efficient Parametric Face Alignment", ICCV 2011

clear; clc; close all
load resultsARSmoothing.mat; %load results

results1 = results(:,2,:,:,:); %Occlusion only
results(:,2,:,:,:) = []; 
results2 = results; %Occlusion + illumination

%%%%%% Just Occlusions
results1 = mean(results1,1);
results1 = squeeze(results1);
affine_ic = results1(:,1)/100;
affine_ECC_ic = results1(:,2)/100;
affine_ic_irls = results1(:,3)/100;
affine_GaborFourier_ic = results1(:,4)/100;
affine_GradientImages_ic = results1(:,5)/100;
affine_GradientCorr_ic = results1(:,6)/100;

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
title('AR: with Smoothing (Occlusion)','Interpreter','tex', 'fontsize', 15)


%%%%%% Occlusions + illumination
results2 = mean(results2,2);
results2 = mean(results2,1);
results2 = squeeze(results2);
affine_ic = results2(:,1)/100;
affine_ECC_ic = results2(:,2)/100;
affine_ic_irls = results2(:,3)/100;
affine_GaborFourier_ic = results2(:,4)/100;
affine_GradientImages_ic = results2(:,5)/100;
affine_GradientCorr_ic = results2(:,6)/100;
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
title('AR: with Smoothing (Occlusion+illumination)','Interpreter','tex', 'fontsize', 15)


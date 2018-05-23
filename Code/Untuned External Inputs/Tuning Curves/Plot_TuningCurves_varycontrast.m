clear all

ROOT     = '/nfs/nhome/live/angus/Documents/SSN_Paper/New_Code/'

X = 'GenerateTuningCurves_varycontrast.mat'; 

%% plot single-cell SI

load([ROOT, X])



figure

subplot(2,2,1)
hold on
set(gca, 'fontsize', 24)
ylabel('Mean Response')
title('Pyramidal Cells')

plot(theta_pE(1:50:end) * 180/pi, RE0_low{1}(1:50:end), 'color', 'r', 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_high{1}(1:50:end), 'color', 'b', 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pE(1:50:end) * 180/pi, RE0_low{2}(1:50:end), 'color', 'r', 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_high{2}(1:50:end), 'color', 'b', 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 3)

plot(theta_pE(1:50:end) * 180/pi, RE0_low{3}(1:50:end), 'color', 'r', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_high{3}(1:50:end), 'color', 'b', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)

legend('Non-spec., low c', 'Non-spec., high c', 'Iso, low c', 'Iso, high c', 'Cross, low c', 'Cross, high c')

plot(theta_pE * 180/pi, RE0_low{1}, 'color', 'r', 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_low{2}, 'color','r', 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_low{3}, 'color', 'r', 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_high{1}, 'color', 'b', 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_high{2}, 'color', 'b', 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_high{3}, 'color', 'b', 'linewidth', 3)



axis([0,360, 0, 0.04])
set(gca, 'xtick', 0:90:360)
set(gca, 'xticklabel', {})
box on

subplot(2,2,2)
hold on
set(gca, 'fontsize', 24)
ylabel('Mean Response')
title('Interneurons')


plot(theta_pI * 180/pi, RI0_low{1}, 'color', 'r', 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_low{1}(1:10:end), 'color', 'r', 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pI * 180/pi, RI0_low{2}, 'color','r', 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_low{2}(1:10:end), 'color', 'r', 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 3)

plot(theta_pI * 180/pi, RI0_low{3}, 'color', 'r', 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_low{3}(1:10:end), 'color', 'r', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)


plot(theta_pI * 180/pi, RI0_high{1}, 'color', 'b', 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_high{1}(1:10:end), 'color', 'b', 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pI * 180/pi, RI0_high{2}, 'color', 'b', 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_high{2}(1:10:end), 'color', 'b', 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 3)

plot(theta_pI * 180/pi, RI0_high{3}, 'color', 'b', 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_high{3}(1:10:end), 'color', 'b', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)


axis([0,360, 0, 0.12])
set(gca, 'xtick', 0:90:360)
set(gca, 'xticklabel', {})
box on

subplot(2,2,3)
hold on
set(gca, 'fontsize', 24)
ylabel('Ratio of Responses')
title('Pyramidal Cells')

plot(theta_pE(1:50:end) * 180/pi, RE0_high{1}(1:50:end) ./ RE0_low{1}(1:50:end), 'color', 'k', 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_high{2}(1:50:end) ./ RE0_low{2}(1:50:end), 'color', 'k', 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi,  RE0_high{3}(1:50:end) ./ RE0_low{3}(1:50:end), 'color', 'k', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)


plot(theta_pE * 180/pi, RE0_high{1} ./ RE0_low{1}, 'color', 'k', 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_high{2} ./  RE0_low{2}, 'color','k', 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_high{3} ./ RE0_low{3}, 'color', 'k', 'linewidth', 3)

legend('Non-spec.', 'Iso', 'Cross')

axis([0,360, 0, 3.5])
set(gca, 'xtick', 0:90:360)
box on

subplot(2,2,4)
hold on
set(gca, 'fontsize', 24)
ylabel('Ratio of responses')
title('Interneurons')

plot(theta_pI * 180/pi, RI0_high{1} ./ RI0_low{1}, 'color', 'k', 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_high{1}(1:10:end) ./ RI0_low{1}(1:10:end), 'color', 'k', 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pI * 180/pi, RI0_high{2} ./  RI0_low{2}, 'color','k', 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_high{2}(1:10:end) ./ RI0_low{2}(1:10:end), 'color', 'k', 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 3)

plot(theta_pI * 180/pi, RI0_high{3} ./ RI0_low{3}, 'color', 'k', 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi,  RI0_high{3}(1:10:end) ./ RI0_low{3}(1:10:end), 'color', 'k', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)

axis([0,360, 0, 3.5])
set(gca, 'xtick', 0:90:360)
box on

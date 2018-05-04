clear all

%ROOT     = '/nfs/nhome/live/angus/Documents/SSN_Paper/New_Code/'

%X = 'GenerateTuningCurves.mat'; 
%X = 'GenerateTuningCurves_matchSI.mat'; 

%% plot single-cell SI

load([ROOT, X])
Xr = rgb2hsv([1,0,0]);
Xb = rgb2hsv([0,0,1]);
Xr(2) = 0.6;
Xb(2) = 0.6;
Xr = hsv2rgb(Xr);
Xb = hsv2rgb(Xb);

Xr2 = rgb2hsv([1,0,0]);
Xb2 = rgb2hsv([0,0,1]);
Xr2(2) = 0.2;
Xb2(2) = 0.2;
Xr2 = hsv2rgb(Xr2);
Xb2 = hsv2rgb(Xb2);



figure

subplot(3,2,1)
hold on
set(gca, 'fontsize', 18)
ylabel('Mean Response')
title('Pyramidal Cells')

plot(theta_pE * 180/pi, RE0_vert{1}, 'color', Xr, 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_vert{1}(1:50:end), 'color', Xr, 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_vert{2}, 'color',Xr2, 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_vert{2}(1:50:end), 'color', Xr2, 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_vert{3}, 'color', 'r', 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_vert{3}(1:50:end), 'color', 'r', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_ang{1}, 'color', Xb, 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_ang{1}(1:50:end), 'color', Xb, 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_ang{2}, 'color', Xb2, 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_ang{2}(1:50:end), 'color', Xb2, 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 2)

plot(theta_pE * 180/pi, RE0_ang{3}, 'color', 'b', 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_ang{3}(1:50:end), 'color', 'b', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)


axis([0,360, 0, 1.1 * max(RE0_vert{3})])
set(gca, 'xtick', 0:90:360)
set(gca, 'xticklabel', {})
box on

subplot(3,2,2)
hold on
set(gca, 'fontsize', 18)
ylabel('Mean Response')
title('Interneurons')


plot(theta_pI * 180/pi, RI0_vert{1}, 'color', Xr, 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_vert{1}(1:10:end), 'color', Xr, 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pI * 180/pi, RI0_vert{2}, 'color',Xr2, 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_vert{2}(1:10:end), 'color', Xr2, 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 3)

plot(theta_pI * 180/pi, RI0_vert{3}, 'color', 'r', 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_vert{3}(1:10:end), 'color', 'r', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)


plot(theta_pI * 180/pi, RI0_ang{1}, 'color', Xb, 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_ang{1}(1:10:end), 'color', Xb, 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pI * 180/pi, RI0_ang{2}, 'color', Xb2, 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_ang{2}(1:10:end), 'color', Xb2, 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 2)

plot(theta_pI * 180/pi, RI0_ang{3}, 'color', 'b', 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_ang{3}(1:10:end), 'color', 'b', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)


axis([0,360, 0, 1.1 * max(RI0_vert{3})])
set(gca, 'xtick', 0:90:360)
set(gca, 'xticklabel', {})
box on

subplot(3,2,3)
hold on
set(gca, 'fontsize', 18)
ylabel('Response std')
title('Pyramidal Cells')


plot(theta_pE * 180/pi, RE0_std_vert{1}, 'color', Xr, 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_std_vert{1}(1:50:end), 'color', Xr, 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_std_vert{2}, 'color',Xr2, 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_std_vert{2}(1:50:end), 'color', Xr2, 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_std_vert{3}, 'color', 'r', 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_std_vert{3}(1:50:end), 'color', 'r', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_std_ang{1}, 'color', Xb, 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_std_ang{1}(1:50:end), 'color', Xb, 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_std_ang{2}, 'color', Xb2, 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_std_ang{2}(1:50:end), 'color', Xb2, 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 2)

plot(theta_pE * 180/pi, RE0_std_ang{3}, 'color', 'b', 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_std_ang{3}(1:50:end), 'color', 'b', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)



axis([0,360, 0, 1.1 * max(RE0_std_vert{3})])
set(gca, 'xtick', 0:90:360)
set(gca, 'xticklabel', {})
box on

subplot(3,2,4)
hold on
set(gca, 'fontsize', 18)
ylabel('Response std')
title('Interneurons')


plot(theta_pI * 180/pi, RI0_std_vert{1}, 'color', Xr, 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_std_vert{1}(1:10:end), 'color', Xr, 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pI * 180/pi, RI0_std_vert{2}, 'color',Xr2, 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_std_vert{2}(1:10:end), 'color', Xr2, 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 3)

plot(theta_pI * 180/pi, RI0_std_vert{3}, 'color', 'r', 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_std_vert{3}(1:10:end), 'color', 'r', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)


plot(theta_pI * 180/pi, RI0_std_ang{1}, 'color', Xb, 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_std_ang{1}(1:10:end), 'color', Xb, 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pI * 180/pi, RI0_std_ang{2}, 'color', Xb2, 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_std_ang{2}(1:10:end), 'color', Xb2, 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 2)

plot(theta_pI * 180/pi, RI0_std_ang{3}, 'color', 'b', 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_std_ang{3}(1:10:end), 'color', 'b', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)


axis([0,360, 0, 1.1 * max(RI0_std_vert{3})])
set(gca, 'xtick', 0:90:360)
set(gca, 'xticklabel', {})
box on

subplot(3,2,5)
hold on
set(gca, 'fontsize', 18)
ylabel('Response Selectivity')
title('Pyramidal Cells')

plot(theta_pE * 180/pi, abs(SI0E{1}), 'color', [0.4, 0.4, 0.4], 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, abs(SI0E{1}(1:50:end)), 'color', [0.4,0.4,0.4], 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pE * 180/pi, abs(SI0E{2}), 'color', [0.8, 0.8, 0.8], 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, abs(SI0E{2}(1:50:end)), 'color', [0.8,0.8,0.8], 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 3)

plot(theta_pE * 180/pi, abs(SI0E{3}), 'k', 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, abs(SI0E{3}(1:50:end)), 'color', 'k', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)

axis([0,360, -0.5, 2.5])
set(gca, 'xtick', 0:90:360)
xlabel('Neuron Angle on Ring')
box on

subplot(3,2,6)
hold on
set(gca, 'fontsize', 18)
ylabel('Response Selectivity')
title('Interneurons')


plot(theta_pI * 180/pi, abs(SI0I{1}), 'color', [0.4, 0.4, 0.4], 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, abs(SI0I{1}(1:10:end)), 'color', [0.4,0.4,0.4], 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pI * 180/pi, abs(SI0I{2}), 'color', [0.8, 0.8, 0.8], 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, abs(SI0I{2}(1:10:end)), 'color', [0.8,0.8,0.8], 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 3)

plot(theta_pI * 180/pi, abs(SI0I{3}), 'k', 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, abs(SI0I{3}(1:10:end)), 'color', 'k', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)


axis([0,360, -0.5, 2.5])
set(gca, 'xtick', 0:90:360)
box on
xlabel('Neuron Angle on Ring') 

%% Plot External and Recurrent Inputs 

figure

for m=1:3

subplot(3,2,2*m-1)
hold on

plot(theta_pE * 180/pi, ExternalE_Input{m}, 'linewidth', 4, 'color', 'k', 'linestyle', '--')
plot(theta_pE * 180/pi, RecurrentEE_Input{m}, 'linewidth', 4, 'color', 'b')
plot(theta_pE * 180/pi, -RecurrentEI_Input{m}, 'linewidth', 4, 'color', 'r')
plot(theta_pE * 180/pi, ExternalE_Input{m} + RecurrentEE_Input{m} - RecurrentEI_Input{m}, 'linewidth', 4, 'color', 'k')

set(gca, 'fontsize', 18)
xlabel('Neuron Angle on Ring')
ylabel('Input Current')
if m==1
legend('External', 'Local E', 'Local I', 'Total', 'Orientation', 'horizontal')
title('Pyramidals (Nonspecific)')
set(gca, 'xticklabel', {})
elseif m==2
title('Pyramidals (Iso-Orientation)')
set(gca, 'xticklabel', {})
elseif m==3
title('Pyramidals (Cross-Orientation)')
end
axis([0,360, -.6,.6])
box on

subplot(3,2, 2*m )
hold on

plot(theta_pI * 180/pi, ExternalI_Input{m}, 'linewidth', 4, 'color', 'k', 'linestyle', '--')
plot(theta_pI * 180/pi, RecurrentIE_Input{m}, 'linewidth', 4, 'color', 'b')
plot(theta_pI * 180/pi, -RecurrentII_Input{m}, 'linewidth', 4, 'color', 'r')
plot(theta_pI * 180/pi, ExternalI_Input{m} + RecurrentIE_Input{m} - RecurrentII_Input{m}, 'linewidth', 4, 'color', 'k')

set(gca, 'fontsize', 18)
xlabel('Neuron Angle on Ring')
ylabel('Input Current')
if m==1
title('Interneurons (Nonspecific)')
set(gca, 'xticklabel', {})
elseif m==2
title('Interneurons (Iso-Orientation)')
set(gca, 'xticklabel', {})
elseif m==3
title('Interneurons (Cross-Orientation)')
end
axis([0,360, -.6,.6])
box on

end

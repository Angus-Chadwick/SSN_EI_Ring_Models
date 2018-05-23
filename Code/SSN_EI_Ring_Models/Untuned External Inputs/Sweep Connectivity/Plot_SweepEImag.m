ROOT     = '/nfs/nhome/live/angus/Documents/SSN_Paper/New_Code/'

X = 'SweepEImag_nonspecificI_squaredSIconvention.mat'; 
Y = 'SweepEImag_isoI_squaredSIconvention.mat'; 
Z = 'SweepEImag_crossI_squaredSIconvention.mat'; 

%% plot single-cell SI

load([ROOT, X])

figure

subplot(2,3,1)


h = imagesc(JIE_mean,JEI_mean, log(sqrt(SI_E')))
set(h, 'AlphaData', SI_E' ~= 0)
set(gca, 'fontsize', 18)
set(gca, 'ydir', 'normal')
set(gca, 'clim', [-3,3])
ylabel('Magnitude of E to I weights (J_{IE})')
xlabel('Magnitude of I to E weights (J_{EI})')
title('Selectivity (Pyramidal Cells)')
h = colorbar;
set(h, 'fontsize', 18)
title(h, 'log SI')

hold on
plot(0.04,0.04, 'linewidth',3, 'linestyle','none', 'marker', '+', 'markersize', 25, 'color', [0.4,0.4,0.4])
plot(0.04,0.04, 'linewidth',3, 'linestyle','none', 'marker', '*', 'markersize', 25, 'color', [0.4,0.4,0.4])
plot(0.04,0.04, 'linewidth',3, 'linestyle','none', 'marker', 'x', 'markersize', 25, 'color', [0.4,0.4,0.4])

subplot(2,3,4)

h = imagesc(JIE_mean,JEI_mean, log(sqrt(SI_I')))
set(h, 'AlphaData', SI_I' ~= 0)
set(gca, 'fontsize', 18)
set(gca, 'ydir', 'normal')
set(gca, 'clim', [-3,3])
ylabel('Magnitude of E to I weights (J_{IE})')
xlabel('Magnitude of I to E weights (J_{EI})')
title('Selectivity (Interneurons)')
h = colorbar;
set(h, 'fontsize', 18)
title(h, 'log SI')

hold on
plot(0.04,0.04, 'linewidth',3, 'linestyle','none', 'marker', '+', 'markersize', 25, 'color', [0.4,0.4,0.4])
plot(0.04,0.04, 'linewidth',3, 'linestyle','none', 'marker', '*', 'markersize', 25, 'color', [0.4,0.4,0.4])
plot(0.04,0.04, 'linewidth',3, 'linestyle','none', 'marker', 'x', 'markersize', 25, 'color', [0.4,0.4,0.4])

load([ROOT, Y])

subplot(2,3,2)


h = imagesc(JIE_mean,JEI_mean, log(sqrt(SI_E')))
set(h, 'AlphaData', SI_E' ~= 0)
set(gca, 'fontsize', 18)
set(gca, 'ydir', 'normal')
set(gca, 'clim', [-3,3])
ylabel('Magnitude of E to I weights (J_{IE})')
xlabel('Magnitude of I to E weights (J_{EI})')
title('Selectivity (Pyramidal Cells)')
h = colorbar;
set(h, 'fontsize', 18)
title(h, 'log SI')


hold on
plot(0.04,0.04, 'linewidth',3, 'linestyle','none', 'marker', '+', 'markersize', 25, 'color', [0.4,0.4,0.4])
plot(0.04,0.04, 'linewidth',3, 'linestyle','none', 'marker', '*', 'markersize', 25, 'color', [0.4,0.4,0.4])
plot(0.04,0.04, 'linewidth',3, 'linestyle','none', 'marker', 'x', 'markersize', 25, 'color', [0.4,0.4,0.4])

subplot(2,3,5)

h = imagesc(JIE_mean,JEI_mean, log(sqrt(SI_I')))
set(h, 'AlphaData', SI_I' ~= 0)
set(gca, 'fontsize', 18)
set(gca, 'ydir', 'normal')
set(gca, 'clim', [-3,3])
ylabel('Magnitude of E to I weights (J_{IE})')
xlabel('Magnitude of I to E weights (J_{EI})')
title('Selectivity (Interneurons)')
h = colorbar;
set(h, 'fontsize', 18)
title(h, 'log SI')

hold on
plot(0.04,0.04, 'linewidth',3, 'linestyle','none', 'marker', '+', 'markersize', 25, 'color', [0.4,0.4,0.4])
plot(0.04,0.04, 'linewidth',3, 'linestyle','none', 'marker', '*', 'markersize', 25, 'color', [0.4,0.4,0.4])
plot(0.04,0.04, 'linewidth',3, 'linestyle','none', 'marker', 'x', 'markersize', 25, 'color', [0.4,0.4,0.4])

load([ROOT, Z])

subplot(2,3,3)


h = imagesc(JIE_mean,JEI_mean, log(sqrt(SI_E')))
set(h, 'AlphaData', SI_E' ~= 0)
set(gca, 'fontsize', 18)
set(gca, 'ydir', 'normal')
set(gca, 'clim', [-3,3])
ylabel('Magnitude of E to I weights (J_{IE})')
xlabel('Magnitude of I to E weights (J_{EI})')
title('Selectivity (Pyramidal Cells)')
h = colorbar;
set(h, 'fontsize', 18)
title(h, 'log SI')

hold on
plot(0.04,0.04, 'linewidth',3, 'linestyle','none', 'marker', '+', 'markersize', 25, 'color', [0.4,0.4,0.4])
plot(0.04,0.04, 'linewidth',3, 'linestyle','none', 'marker', '*', 'markersize', 25, 'color', [0.4,0.4,0.4])
plot(0.04,0.04, 'linewidth',3, 'linestyle','none', 'marker', 'x', 'markersize', 25, 'color', [0.4,0.4,0.4])

subplot(2,3,6)

h = imagesc(JIE_mean,JEI_mean, log(sqrt(SI_I')))
set(h, 'AlphaData', SI_I' ~= 0)
set(gca, 'fontsize', 18)
set(gca, 'ydir', 'normal')
set(gca, 'clim', [-3,3])
ylabel('Magnitude of E to I weights (J_{IE})')
xlabel('Magnitude of I to E weights (J_{EI})')
title('Selectivity (Interneurons)')
h = colorbar;
set(h, 'fontsize', 18)
title(h, 'log SI')

hold on
plot(0.04,0.04, 'linewidth',3, 'linestyle','none', 'marker', '+', 'markersize', 25, 'color', [0.4,0.4,0.4])
plot(0.04,0.04, 'linewidth',3, 'linestyle','none', 'marker', '*', 'markersize', 25, 'color', [0.4,0.4,0.4])
plot(0.04,0.04, 'linewidth',3, 'linestyle','none', 'marker', 'x', 'markersize', 25, 'color', [0.4,0.4,0.4])

colormap copper

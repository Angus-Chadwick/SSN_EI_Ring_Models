clear all

%ROOT     = '/nfs/nhome/live/angus/Documents/SSN_Paper/New_Code/'
ROOT     = '/mnt/data/angusc/Data/Tuned External Inputs/Sweep Connectivity/';

X = 'SweepEIwidth_TunedExtIn.mat'; 

%% plot single-cell SI

load([ROOT, X])

figure

axes('Position', [0.075, 0.1, 0.4, 0.35])

h = imagesc(([1:Nvals] - 26) /20, ([1:Nvals] - 26) /20, log(sqrt(SI_E')))
set(h, 'AlphaData', SI_E' ~= 0)
set(gca, 'fontsize', 18)
set(gca, 'ydir', 'normal')
set(gca, 'clim', [-3,3])
ylabel('Tuning of E to I weights (k_{IE})')
xlabel('Tuning of I to E weights (k_{EI})')
title('Selectivity (Pyramidal Cells)')
axis([(1-26)/20, (Nvals - 26)/20, 0, (Nvals - 26)/20])
h = colorbar;
set(h, 'fontsize', 18)
title(h, 'log SI')

hold on
plot(0.5,0.5, 'linewidth',3, 'linestyle','none', 'marker', '+', 'markersize', 25, 'color', [0.418,0.418,0.418])
plot(0.0,0.0, 'linewidth',3, 'linestyle','none', 'marker', '*', 'markersize', 25, 'color', [0.418,0.418,0.418])
plot(-0.5,0.5, 'linewidth',3, 'linestyle','none', 'marker', 'x', 'markersize', 25, 'color', [0.418,0.418,0.418])

axes('Position', [0.575, 0.1, 0.4, 0.35])

h = imagesc(([1:Nvals] - 26) /20, ([1:Nvals] - 26) /20, log(sqrt(SI_I')))
set(h, 'AlphaData', SI_I' ~= 0)
set(gca, 'fontsize', 18)
set(gca, 'ydir', 'normal')
set(gca, 'clim', [-3,3])
ylabel('Tuning of E to I weights (k_{IE})')
xlabel('Tuning of I to E weights (k_{EI})')
title('Selectivity (Interneurons)')
axis([(1-26)/20, (Nvals - 26)/20, 0, (Nvals - 26)/20])
h = colorbar;
set(h, 'fontsize', 18)
title(h, 'log SI')

hold on
plot(0.5,0.5, 'linewidth',3, 'linestyle','none', 'marker', '+', 'markersize', 25, 'color', [0.418,0.418,0.418])
plot(0.0,0.0, 'linewidth',3, 'linestyle','none', 'marker', '*', 'markersize', 25, 'color', [0.418,0.418,0.418])
plot(-0.5,0.5, 'linewidth',3, 'linestyle','none', 'marker', 'x', 'markersize', 25, 'color', [0.418,0.418,0.418])

%colormap copper
ROOT     = '/nfs/nhome/live/angus/Documents/SSN_Paper/New_Code/'

X = 'SweepEIwidth_squaredSIconvention.mat'; 
%X = 'SweepEIwidth_squaredSIconvention_TunedExtIn.mat'; 

%% plot single-cell SI

load([ROOT, X])

figure

axes('Position', [0.05, 0.1, 0.45, 0.3])

h = imagesc(([1:Nvals] - 26) /20, ([1:Nvals] - 26) /20, sqrt(SI_E'))
set(h, 'AlphaData', SI_E' ~= 0)
set(gca, 'fontsize', 24)
set(gca, 'ydir', 'normal')
set(gca, 'clim', [0,8])
ylabel('Tuning of E to I weights (k_{IE})')
xlabel('Tuning of I to E weights (k_{EI})')
title('Selectivity (Pyramidal Cells)')
axis([(1-26)/20, (Nvals - 26)/20, 0, (Nvals - 26)/20])
h = colorbar;
set(h, 'fontsize', 24)
title(h, 'SI')

hold on
plot(0.3,0.3, 'linewidth',3, 'linestyle','none', 'marker', '+', 'markersize', 25, 'color', [0.418,0.418,0.418])
plot(0.0,0.0, 'linewidth',3, 'linestyle','none', 'marker', '*', 'markersize', 25, 'color', [0.418,0.418,0.418])
plot(-0.3,0.3, 'linewidth',3, 'linestyle','none', 'marker', 'x', 'markersize', 25, 'color', [0.418,0.418,0.418])

axes('Position', [0.55, 0.1, 0.45, 0.3])

h = imagesc(([1:Nvals] - 26) /20, ([1:Nvals] - 26) /20, sqrt(SI_I'))
set(h, 'AlphaData', SI_I' ~= 0)
set(gca, 'fontsize', 24)
set(gca, 'ydir', 'normal')
set(gca, 'clim', [0,15])
ylabel('Tuning of E to I weights (k_{IE})')
xlabel('Tuning of I to E weights (k_{EI})')
title('Selectivity (Interneurons)')
axis([(1-26)/20, (Nvals - 26)/20, 0, (Nvals - 26)/20])
h = colorbar;
set(h, 'fontsize', 24)
title(h, 'SI')

hold on
plot(0.3,0.3, 'linewidth',3, 'linestyle','none', 'marker', '+', 'markersize', 25, 'color', [0.418,0.418,0.418])
plot(0.0,0.0, 'linewidth',3, 'linestyle','none', 'marker', '*', 'markersize', 25, 'color', [0.418,0.418,0.418])
plot(-0.3,0.3, 'linewidth',3, 'linestyle','none', 'marker', 'x', 'markersize', 25, 'color', [0.418,0.418,0.418])

colormap copper
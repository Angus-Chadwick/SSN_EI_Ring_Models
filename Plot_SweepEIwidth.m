ROOT     = '/nfs/nhome/live/angus/Documents/SSN_Paper/New_Code/'

X = 'SweepEIwidth_squaredSIconvention.mat'; 

%% plot single-cell SI

load([ROOT, X])

figure

subplot(1,2,1)

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

subplot(1,2,2)

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

%% plot pop SI

figure 
h = imagesc(([1:Nvals] - 26) /20, ([1:Nvals] - 26) /20, sqrt(SItot_E') .* (SItot_E' > 0))
set(h, 'AlphaData', ~isnan(SItot_E'))
set(gca, 'fontsize', 24)
set(gca, 'ydir', 'normal')
set(gca, 'clim', [0,8])
ylabel('Tuning of E to I weights (k_{IE})')
xlabel('Tuning of I to E weights (k_{EI})')
title('Selectivity (Pyramidal Cell Population)')
axis([(1-26)/20, (Nvals - 26)/20, 0, (Nvals - 26)/20])
h = colorbar;
set(h, 'fontsize', 24)
title(h, 'SI')


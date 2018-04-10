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

%% Plot population SI


load([ROOT, X])

figure

subplot(1,3,1)

h = imagesc(JEI_mean,JIE_mean, log(sqrt(SItot_E .* (SItot_E > 0))))
set(h, 'AlphaData', ~isnan(SItot_E))
set(gca, 'fontsize', 18)
set(gca, 'ydir', 'normal')
xlabel('E to I coupling strength (J_{IE})')
ylabel('I to E coupling strength (J_{EI})')
title('Selectivity (Pyramidal Population)')
h = colorbar;
set(h, 'fontsize', 18)
title(h, 'log SI')

load([ROOT, Y])



subplot(1,3,2)

h = imagesc(JEI_mean,JIE_mean, log(sqrt(SItot_E .* (SItot_E > 0))))
set(h, 'AlphaData', ~isnan(SItot_E))
set(gca, 'fontsize', 18)
set(gca, 'ydir', 'normal')
xlabel('E to I coupling strength (J_{IE})')
ylabel('I to E coupling strength (J_{EI})')
title('Selectivity (Pyramidal Population)')
h = colorbar;
set(h, 'fontsize', 18)
title(h, 'log SI')

load([ROOT, Z])



subplot(1,3,3)

h = imagesc(JEI_mean,JIE_mean, log(sqrt(SItot_E .* (SItot_E > 0))))
set(h, 'AlphaData', ~isnan(SItot_E))
set(gca, 'fontsize', 18)
set(gca, 'ydir', 'normal')
xlabel('E to I coupling strength (J_{IE})')
ylabel('I to E coupling strength (J_{EI})')
title('Selectivity (Pyramidal Population)')
h = colorbar;
set(h, 'fontsize', 18)
title(h, 'log SI')
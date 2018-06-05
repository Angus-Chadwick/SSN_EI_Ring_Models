% Plot Activity, covariance and spectra


%% Placement Parameters

UL = 0.03;  % upper left figure

boxwidth = 0.2;  % total height of a panel
boxheight = 0.2;  % total width of a panel
boxspace_v = 0.05; % vertical spacing between panels

boxspace_h = 0.05;

inset_size = 0.05;

%% Plot Mean-Eigenvector Overlaps

% Nonspecific-Inhibition

ROOT     = '/mnt/data/angusc/Data/Tuned External Inputs/Covariance/';

X =        'CovarianceData_NSI_TunedExtIn.mat';
load([ROOT, X])

RE_covtot1 = RE_covtot1(1:1000,1:1000);
RE_covtot2 = RE_covtot2(1:1000,1:1000);

[V,D] = eig((0.5 * (RE_covtot1(1:1000,1:1000) + RE_covtot2(1:1000,1:1000) )));

axes('Position', [UL + 1.9 * boxwidth, 1-UL - boxheight, 0.6 * boxwidth, boxheight])

hold on

E_Spect_Overlap = flipud(abs(V' * diff(RE0')' ./ (norm(diff(RE0')))));

plot(E_Spect_Overlap(1), 'linewidth', 4, 'marker', 'x', 'markersize', 8, 'linestyle', 'none', 'color', 'r')
plot(2, E_Spect_Overlap(2), 'linewidth', 4, 'marker', 'x', 'markersize', 8, 'linestyle', 'none', 'color', 'b')
plot(3:1000, E_Spect_Overlap(3:end), 'linewidth', 4, 'marker', 'x', 'markersize', 8, 'linestyle', 'none', 'color', 'k')

hold on
plot([-10,1010], [1,1], 'color', 'k', 'linewidth', 2, 'linestyle', '--')
plot([-10,1010], [0,0], 'color', 'k', 'linewidth', 2, 'linestyle', '--')

set(gca,'fontsize', 10)
title('Nonspecific Inhibition')
xlabel('Principal Component')
ylabel('Overlap with Signal Vector')
axis([0.5,10.5,-0.1,1.1])
box on

% Iso-Inhibition

ROOT     = '/mnt/data/angusc/Data/Tuned External Inputs/Covariance/';

X =        'CovarianceData_IsoI_TunedExtIn_Batch1.mat';
load([ROOT, X])

X =        'CovarianceData_IsoI_TunedExtIn_PooledCovsandmeans.mat';
load([ROOT, X])

RE_covtot1 = RE_covtot1(1:1000,1:1000);
RE_covtot2 = RE_covtot2(1:1000,1:1000);

[V,D] = eig((0.5 * (RE_covtot1(1:1000,1:1000) + RE_covtot2(1:1000,1:1000) )));

axes('Position', [UL + 2.8 * boxwidth, 1-UL - boxheight, 0.6 * boxwidth, boxheight])

hold on

E_Spect_Overlap = flipud(abs(V' * diff(RE0')' ./ (norm(diff(RE0')))));

plot(E_Spect_Overlap(1), 'linewidth', 4, 'marker', 'x', 'markersize', 8, 'linestyle', 'none', 'color', 'r')
plot(2, E_Spect_Overlap(2), 'linewidth', 4, 'marker', 'x', 'markersize', 8, 'linestyle', 'none', 'color', 'b')
plot(3:1000, E_Spect_Overlap(3:end), 'linewidth', 4, 'marker', 'x', 'markersize', 8, 'linestyle', 'none', 'color', 'k')

hold on
plot([-10,1010], [1,1], 'color', 'k', 'linewidth', 2, 'linestyle', '--')
plot([-10,1010], [0,0], 'color', 'k', 'linewidth', 2, 'linestyle', '--')

set(gca,'fontsize', 10)
title('Iso-Orienation Inhibition')
xlabel('Principal Component')
ylabel('Overlap with Signal Vector')
axis([0.5,10.5,-0.1,1.1])
box on

% Cross-Inhibition

ROOT     = '/mnt/data/angusc/Data/Tuned External Inputs/Covariance/';

X =        'CovarianceData_CrossI_TunedExtIn.mat';
load([ROOT, X])

RE_covtot1 = RE_covtot1(1:1000,1:1000);
RE_covtot2 = RE_covtot2(1:1000,1:1000);

[V,D] = eig((0.5 * (RE_covtot1(1:1000,1:1000) + RE_covtot2(1:1000,1:1000) )));

axes('Position', [UL + 3.7 * boxwidth, 1-UL - boxheight, 0.6 * boxwidth, boxheight])

hold on

E_Spect_Overlap = flipud(abs(V' * diff(RE0')' ./ (norm(diff(RE0')))));

plot(E_Spect_Overlap(1), 'linewidth', 4, 'marker', 'x', 'markersize', 8, 'linestyle', 'none', 'color', 'r')
plot(2, E_Spect_Overlap(2), 'linewidth', 4, 'marker', 'x', 'markersize', 8, 'linestyle', 'none', 'color', 'b')
plot(3:1000, E_Spect_Overlap(3:end), 'linewidth', 4, 'marker', 'x', 'markersize', 8, 'linestyle', 'none', 'color', 'k')

hold on
plot([-10,1010], [1,1], 'color', 'k', 'linewidth', 2, 'linestyle', '--')
plot([-10,1010], [0,0], 'color', 'k', 'linewidth', 2, 'linestyle', '--')

set(gca,'fontsize', 10)
title('Cross-Orientation Inhibition')
xlabel('Principal Component')
ylabel('Overlap with Signal Vector')
axis([0.5,10.5,-0.1,1.1])
box on

%% Plot Selectivity Sweep and Decompositino

ROOT     = '/mnt/data/angusc/Data/Tuned External Inputs/Sweep Connectivity/';
X = 'SweepEIwidth_TunedExtIn.mat'; 
load([ROOT, X])

axes('Position', [UL + 0.2 * boxwidth, 1-UL - 2.4 * boxheight, 1.3 * boxwidth, 1.0 * boxheight])

h = imagesc(([1:Nvals] - 26) /20, ([1:Nvals] - 26) /20, log(sqrt(SItot_E_epsridge') .* (SItot_E_epsridge' > 0)))
set(h, 'AlphaData', ~isnan(SItot_E_epsridge'))
set(gca, 'fontsize', 10)
set(gca, 'ydir', 'normal')
set(gca, 'clim', [1,4])
ylabel('Tuning of E to I weights (k_{IE})')
xlabel('Tuning of I to E weights (k_{EI})')
title('Decoder Performance')
axis([(1-26)/20, (Nvals - 26)/20, 0, (Nvals - 26)/20])
h = colorbar;
title(h, 'log SI')
set(h, 'fontsize', 10)
box on


ROOT     = '/mnt/data/angusc/Data/Tuned External Inputs/Covariance/';

X = 'DprimeFactors_VaryNonspecificItoIsoI.mat';
load([ROOT, X])

axes('Position', [UL + 1.9 * boxwidth, 1-UL - 2.4 * boxheight, 0.6 * boxwidth, 1.0 * boxheight])

hold on
plot(-1:0.1:0, fliplr(sqrt(NrmSq)), 'linewidth', 4, 'color', [0.3,0.3,0.3])
set(gca, 'fontsize', 10)
set(gca, 'xtick', [-1,0,1]);
set(gca, 'xticklabel', {'Iso', 'N.S.', 'Ortho'})
xlabel('Network Type')
ylabel('Signal (||\Delta r||)')

X = 'DprimeFactors_VaryNonspecificItoCrossI.mat';
load([ROOT, X])

plot(0:0.1:1, sqrt(NrmSq), 'linewidth', 4, 'color', [0.7,0.7,0.7])
box on
axis([-1,1, -0.1,0.5])

X = 'DprimeFactors_VaryNonspecificItoIsoI.mat';
load([ROOT, X])

axes('Position', [UL + 2.8 * boxwidth, 1-UL - 2.4 * boxheight, 0.6 * boxwidth, 1.0 * boxheight])

hold on
plot(-1:0.1:0, fliplr(sqrt(TrD)), 'linewidth', 4, 'color', [0.3,0.3,0.3])
set(gca, 'fontsize', 10)
set(gca, 'xtick', [-1,0,1]);
set(gca, 'xticklabel', {'Iso', 'N.S.', 'Ortho'})
xlabel('Network Type')
ylabel('Noise (\surd{ Trace \Sigma })')
box on

X = 'DprimeFactors_VaryNonspecificItoCrossI.mat';
load([ROOT, X])

plot(0:0.1:1, sqrt(TrD), 'linewidth', 4, 'color', [0.7,0.7,0.7])
axis([-1,1, -0.1,0.5])


axes('Position', [UL + 3.7 * boxwidth, 1-UL - 2.4 * boxheight, 0.6 * boxwidth, 1.0 * boxheight])

X = 'DprimeFactors_VaryNonspecificItoIsoI.mat';
load([ROOT, X])

hold on
plot(-1:0.1:0, fliplr(1./sqrt(SNRalign)), 'linewidth', 4, 'color', [0.3,0.3,0.3])

X = 'DprimeFactors_VaryNonspecificItoCrossI.mat';
load([ROOT, X])

plot(0:0.1:1, 1./sqrt(SNRalign), 'linewidth', 4, 'color', [0.7,0.7,0.7])
set(gca, 'fontsize', 10)
set(gca, 'xtick', [-1,0,1]);
set(gca, 'xticklabel', {'Iso', 'N.S.', 'Ortho'})
xlabel('Network Type')
ylabel('Signal-Noise Alignment')
box on
axis([-1,1, 0,0.15])


%% Plot Sweeps over EI mags

ROOT     = '/mnt/data/angusc/Data/Tuned External Inputs/Sweep Connectivity/';
         
X1 = 'SweepEImag_NonspecificI_TunedExtIn.mat'; 
X2 = 'SweepEImag_IsoI_TunedExtIn.mat'; 
X3 = 'SweepEImag_CrossI_TunedExtIn.mat'; 

load([ROOT, X1])

JEI_mean = JEI_mean(:,1);
JIE_mean = JIE_mean(1,:);

axes('Position', [UL + 0.3 * boxwidth, 1-UL - 4.3 * boxheight, 1 * boxwidth, 1.2 * boxheight])

h = imagesc(JEI_mean,JIE_mean, log(sqrt(SItot_E_epsridge' .* (SItot_E_epsridge' > 0))))
set(h, 'AlphaData', ~isnan(SItot_E_epsridge'))
set(gca, 'fontsize', 10)
set(gca, 'clim', [1,4])
set(gca, 'ydir', 'normal')
ylabel('E to I coupling strength (J_{IE})')
xlabel('I to E coupling strength (J_{EI})')
title('Decoder Performance')
h = colorbar;
set(h, 'fontsize', 10)
title(h, 'log SI')

load([ROOT, X2])

JEI_mean = JEI_mean(:,1);
JIE_mean = JIE_mean(1,:);

axes('Position', [UL + 1.8  * boxwidth, 1-UL - 4.3 * boxheight, 1 * boxwidth, 1.2 * boxheight])

h = imagesc(JEI_mean,JIE_mean, log(sqrt(SItot_E_epsridge' .* (SItot_E_epsridge' > 0))))
set(h, 'AlphaData', ~isnan(SItot_E_epsridge))
set(gca, 'fontsize', 10)
set(gca, 'clim', [1,4])
set(gca, 'ydir', 'normal')
ylabel('E to I coupling strength (J_{IE})')
xlabel('I to E coupling strength (J_{EI})')
title('Decoder Performance')
h = colorbar;
set(h, 'fontsize', 10)
title(h, 'log SI')


load([ROOT, X3])

JEI_mean = JEI_mean(:,1);
JIE_mean = JIE_mean(1,:);

axes('Position', [UL + 3.3 * boxwidth, 1-UL - 4.3 * boxheight, 1 * boxwidth, 1.2 * boxheight])

h = imagesc(JEI_mean,JIE_mean, log(sqrt(SItot_E_epsridge' .* (SItot_E_epsridge' > 0))))
set(h, 'AlphaData', ~isnan(SItot_E_epsridge'))
set(gca, 'fontsize', 10)
set(gca, 'clim', [1,4])
set(gca, 'ydir', 'normal')
ylabel('E to I coupling strength (J_{IE})')
xlabel('I to E coupling strength (J_{EI})')
title('Decoder Performance')
h = colorbar;
set(h, 'fontsize', 10)
title(h, 'log SI')

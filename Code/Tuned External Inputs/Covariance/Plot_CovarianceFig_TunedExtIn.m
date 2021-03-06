% Plot Activity, covariance and spectra

%% Load Data

clear all

ROOT     = '/mnt/data/angusc/Data/Tuned External Inputs/Covariance/';

X =        'CovarianceData_IsoI_TunedExtIn_Batch1.mat';
load([ROOT, X])

X =        'CovarianceData_IsoI_TunedExtIn_PooledCovsandmeans.mat'; % REDO WITH FULL RE0 TOO
load([ROOT, X])

RE_covtot1 = RE_covtot1(1:1000,1:1000);
RE_covtot2 = RE_covtot2(1:1000,1:1000);

[V,D] = eig((0.5 * (RE_covtot1(1:1000,1:1000) + RE_covtot2(1:1000,1:1000) )));

%% Placement Parameters

UL = 0.03;  % upper left figure

boxwidth = 0.2;  % total height of a panel
boxheight = 0.2;  % total width of a panel
boxspace_v = 0.05; % vertical spacing between panels

boxspace_h = 0.05;

inset_size = 0.05;

%% Plot Raw Data Samples

axes('Position', [UL, 1-UL - boxheight/2 + boxspace_v/2, 0.8 * boxwidth, boxheight/2 - boxspace_v/2])

imagesc((1:Nt) / 1000, 1:1000,rE); 
h = colorbar;

cmax = 0.8 * max(max(rE));
cmin = 0;
set(gca, 'clim', [cmin, cmax])

set(gca, 'fontsize', 10)
set(gca, 'xtick', [])
set(gca, 'ytick', [])
set(h, 'fontsize', 10)
ylabel('E cells')
title('Population Activity')
box on

axes('Position', [UL, 1-UL - boxheight, 0.8 * boxwidth, boxheight/2 - boxspace_v/2])

imagesc((1:Nt)/1000, 1:200, rI)  % set size of this window to be small
h = colorbar;

cmax = 0.8 * max(max(rI));
cmin = 0;
set(gca, 'clim', [cmin, cmax])

set(gca, 'fontsize', 10)
set(gca, 'xtick', [])
set(gca, 'ytick', [])
set(h, 'fontsize', 10)
xlabel('Time (s)')
ylabel('I cells')
box on

%% Plot Covariance Matrix

axes('Position', [UL + 0.8  * boxwidth + 0.5 * boxspace_h,1-UL - boxheight, 0.8 * boxwidth, boxheight])

CovMat = RE_covtot2 - diag(diag(RE_covtot2));

imagesc(CovMat)
h = colorbar;

hold on

plot(1000:1000, 0:1200, 'color', [1,1,1])
plot(0:1200, 1000:1000, 'color', [1,1,1])

cmax = max(max(CovMat(1:1000,1:1000)));
cmin = -cmax;


set(gca, 'fontsize', 10)
set(gca, 'clim', [cmin, cmax])
set(h, 'fontsize', 10)
set(gca, 'xtick', [])
set(gca, 'ytick', [])
xlabel('Neuron')
ylabel('Neuron')
title('Covariance')
%colormap copper
box on

%% Plot Eigenspectrum

axes('Position', [UL + 1.9 * boxwidth, 1-UL - boxheight, 0.6 * boxwidth, boxheight])
    
hold on

E_Spect = flipud(diag(D)) / trace(D) * 100;

plot(1,E_Spect(1), 'color', 'r', 'marker', 'x', 'markersize', 8, 'linestyle', 'none', 'linewidth', 4)
plot(2,E_Spect(2), 'color', 'b', 'marker', 'x', 'markersize', 8, 'linestyle', 'none', 'linewidth', 4)
plot(3:1000,E_Spect(3:end), 'color', 'k', 'marker', 'x', 'markersize', 8, 'linestyle', 'none', 'linewidth', 4)

axis([0.5,10.5, 0,2.5])
set(gca, 'fontsize', 10)
xlabel('Principal Component')
ylabel('Variance Explained (%)')
box on

axes('Position',[UL + 1.9 * boxwidth + 0.6 * boxwidth - 1.25 * inset_size,1-UL - boxheight + boxheight - 1.25 * inset_size, inset_size, inset_size])
plot(flipud(diag(D)) / trace(D) * 100, 'color', 'k', 'linewidth', 4)
set(gca, 'fontsize', 10)
xlabel('Princ. Comp.')
ylabel('Var. Expl. (%)')
axis([-10, 1010, -0.5, 2.5])
box on

axes('Position', [UL + 2.8 * boxwidth, 1-UL - boxheight/2 + boxspace_v/2, 0.6 * boxwidth, boxheight/2 - boxspace_v/2])

plot([1:50:1000,1000],V([1:50:1000,1000],end), 'linewidth', 4, 'marker', 'x', 'markersize', 8, 'color', 'r', 'linestyle', 'none')
axis([0,1000, -0.125,0.125])
set(gca, 'fontsize', 10)
set(gca, 'xticklabel', {})
box on

legend({'1st P.C.'}, 'Location', 'SouthEast')
legend boxoff

axes('Position', [UL + 2.8 * boxwidth, 1-UL - boxheight, 0.6 * boxwidth, boxheight/2 - boxspace_v/2])

plot([1:50:1000,1000],-V([1:50:1000,1000],end-1), 'linewidth', 4, 'marker', 'x', 'markersize', 8, 'color', 'b', 'linestyle', 'none')
axis([0,1000, -0.125,0.125])
set(gca, 'fontsize', 10)
xlabel('Pyramidal cell')
ylabel('Loading Factor')
box on

legend({'2nd P.C.'}, 'Location', 'SouthEast')
legend boxoff

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
xlabel('Principal Component')
ylabel('Overlap with Mean Vector')
axis([0.5,10.5,-0.1,1.1])
box on

%% Plot Selectivity Sweep and Decompositino

ROOT     = '/mnt/data/angusc/Data/Tuned External Inputs/Sweep Connectivity/';
X = 'SweepEIwidth_TunedExtIn.mat'; 
load([ROOT, X])

axes('Position', [UL + 0.1 * boxwidth, 1-UL - 2.5 * boxheight, 1.5 * boxwidth, 1.0 * boxheight])

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
plot(-1:0.1:0, fliplr(NrmSq), 'linewidth', 4, 'color', [0.3,0.3,0.3])
set(gca, 'fontsize', 10)
set(gca, 'xtick', [-1,0,1]);
set(gca, 'xticklabel', {'Iso', 'N.S.', 'Ortho'})
xlabel('Network Type')
ylabel('Signal (||\Delta r||^2)')

X = 'DprimeFactors_VaryNonspecificItoCrossI.mat';
load([ROOT, X])

plot(0:0.1:1, NrmSq, 'linewidth', 4, 'color', [0.7,0.7,0.7])
box on
axis([-1,1, -0.1,0.5])

X = 'DprimeFactors_VaryNonspecificItoIsoI.mat';
load([ROOT, X])

axes('Position', [UL + 2.8 * boxwidth, 1-UL - 2.4 * boxheight, 0.6 * boxwidth, 1.0 * boxheight])

hold on
plot(-1:0.1:0, fliplr(TrD), 'linewidth', 4, 'color', [0.3,0.3,0.3])
set(gca, 'fontsize', 10)
set(gca, 'xtick', [-1,0,1]);
set(gca, 'xticklabel', {'Iso', 'N.S.', 'Ortho'})
xlabel('Network Type')
ylabel('Noise (Trace \Sigma)')
box on

X = 'DprimeFactors_VaryNonspecificItoCrossI.mat';
load([ROOT, X])

plot(0:0.1:1, TrD, 'linewidth', 4, 'color', [0.7,0.7,0.7])
axis([-1,1, -0.1,0.5])


axes('Position', [UL + 3.7 * boxwidth, 1-UL - 2.4 * boxheight, 0.6 * boxwidth, 1.0 * boxheight])

X = 'DprimeFactors_VaryNonspecificItoIsoI.mat';
load([ROOT, X])

hold on
plot(-1:0.1:0, fliplr(1./SNRalign), 'linewidth', 4, 'color', [0.3,0.3,0.3])

X = 'DprimeFactors_VaryNonspecificItoCrossI.mat';
load([ROOT, X])

plot(0:0.1:1, 1./SNRalign, 'linewidth', 4, 'color', [0.7,0.7,0.7])
set(gca, 'fontsize', 10)
set(gca, 'xtick', [-1,0,1]);
set(gca, 'xticklabel', {'Iso', 'N.S.', 'Ortho'})
xlabel('Network Type')
ylabel('Signal-Noise Alignment')
box on
axis([-1,1, 0,0.025])


%% Plot Sweeps over EI mags

ROOT     = '/mnt/data/angusc/Data/Tuned External Inputs/Sweep Connectivity/';
         
X1 = 'SweepEImag_NonspecificI_TunedExtIn.mat'; 
X2 = 'SweepEImag_IsoI_TunedExtIn.mat'; 
X3 = 'SweepEImag_CrossI_TunedExtIn.mat'; 

load([ROOT, X1])

JEI_mean = JEI_mean(:,1);
JIE_mean = JIE_mean(1,:);

axes('Position', [UL + 0.2 * boxwidth, 1-UL - 4.3 * boxheight, 1.2 * boxwidth, 1.2 * boxheight])

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

axes('Position', [UL + 1.7  * boxwidth, 1-UL - 4.3 * boxheight, 1.2 * boxwidth, 1.2 * boxheight])

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

axes('Position', [UL + 3.2 * boxwidth, 1-UL - 4.3 * boxheight, 1.2 * boxwidth, 1.2 * boxheight])

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


%% Save Figure 

print('CovarianceFig', '-dsvg')
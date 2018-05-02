% Plot Activity, covariance and spectra

%% Load Data

clear all

ROOT     = '/nfs/nhome/live/angus/Documents/SSN_Paper/New_Code/'

X1 = 'CovarianceFig_Data.mat'; 
X2 = 'CovarianceFig_Data_2ndbatch.mat'; 

load([ROOT, X1])
Rtot0_cov1 = Rtot0_cov;
RE01 = RE0;

load([ROOT, X2])
Rtot0_cov2 = Rtot0_cov;
RE02 = RE0;

RE0 = (RE01 + RE02)/2;
Rtot0_cov = (Rtot0_cov1 + Rtot0_cov2)/2;

RE_covtot1 = Rtot0_cov(1:1000,1:1000,1);
RE_covtot2 = Rtot0_cov(1:1000,1:1000,2);

[V,D] = eig((0.5 * (RE_covtot1 + RE_covtot2 )));
%% Placement Parameters

UL = 0.05;  % upper left figure

boxwidth = 0.20;  % total height of a panel
boxheight = 0.20;  % total width of a panel
boxspace_v = 0.01; % vertical spacing between panels

boxspace_h = 0.00;

inset_size = 0.05;

%% Plot Raw Data Samples

axes('Position', [UL, 1-UL - boxheight/2 + boxspace_v/2, boxwidth, boxheight/2 - boxspace_v/2])

imagesc((1:Nt) / 1000, 1:1000,rE); 
h = colorbar;

cmax = 0.8 * max(max(rE));
cmin = 0;
set(gca, 'clim', [cmin, cmax])

set(gca, 'fontsize', 24)
set(gca, 'xtick', [])
set(gca, 'ytick', [])
set(h, 'fontsize', 24)
ylabel('E cells')
title('Population Activity')
box on

axes('Position', [UL, 1-UL - boxheight, boxwidth, boxheight/2 - boxspace_v/2])

imagesc((1:Nt)/1000, 1:200, rI)  % set size of this window to be small
h = colorbar;

cmax = 0.8 * max(max(rI));
cmin = 0;
set(gca, 'clim', [cmin, cmax])

set(gca, 'fontsize', 24)
set(gca, 'xtick', [])
set(gca, 'ytick', [])
set(h, 'fontsize', 24)
xlabel('Time (s)')
ylabel('I cells')
box on

%% Plot Covariance Matrix

axes('Position', [UL + 0.9 * boxwidth + boxspace_h,1-UL - boxheight, boxwidth, boxheight])

CovMat = Rtot0_cov(:,:,2) - diag(diag(Rtot0_cov(:,:,2)));

imagesc(CovMat)
h = colorbar;

hold on

plot(1000:1000, 0:1200, 'color', [1,1,1])
plot(0:1200, 1000:1000, 'color', [1,1,1])

cmax = max(max(CovMat(1:1000,1:1000)));
cmin = -cmax;


set(gca, 'fontsize', 24)
set(gca, 'clim', [cmin, cmax])
set(h, 'fontsize', 24)
set(gca, 'xtick', [])
set(gca, 'ytick', [])
xlabel('Neuron')
ylabel('Neuron')
title('Covariance')
colormap copper
box on

%% Plot Eigenspectrum

axes('Position', [UL + 1.9 * boxwidth + 2 * boxspace_h, 1-UL - boxheight, 0.6 * boxwidth, boxheight])
    
hold on

E_Spect = flipud(diag(D)) / trace(D) * 100;

plot(1,E_Spect(1), 'color', 'r', 'marker', 'x', 'markersize', 10, 'linestyle', 'none', 'linewidth', 4)
plot(2,E_Spect(2), 'color', 'b', 'marker', 'x', 'markersize', 10, 'linestyle', 'none', 'linewidth', 4)
plot(3:1000,E_Spect(3:end), 'color', 'k', 'marker', 'x', 'markersize', 10, 'linestyle', 'none', 'linewidth', 4)

axis([0.5,10.5, 0,2.5])
set(gca, 'fontsize', 24)
xlabel('Principal Component')
ylabel('Variance Explained (%)')
box on

axes('Position',[UL + 1.9 * boxwidth + 2 * boxspace_h + 0.6 * boxwidth - 1.25 * inset_size,1-UL - boxheight + boxheight - 1.25 * inset_size, inset_size, inset_size])
plot(flipud(diag(D)) / trace(D) * 100, 'color', 'k', 'linewidth', 4)
set(gca, 'fontsize', 18)
xlabel('Princ. Comp.')
ylabel('Var. Expl. (%)')
axis([-10, 1010, -0.5, 2.5])
box on

axes('Position', [UL + 2.8 * boxwidth + 5 * boxspace_h, 1-UL - boxheight/2 + boxspace_v/2, 0.6 * boxwidth, boxheight/2 - boxspace_v/2])

plot([1:50:1000,1000],-V([1:50:1000,1000],end), 'linewidth', 4, 'marker', 'x', 'markersize', 10, 'color', 'r', 'linestyle', 'none')
axis([0,1000, -0.125,0.125])
set(gca, 'fontsize', 24)
set(gca, 'xticklabel', {})
box on

legend({'1st P.C.'}, 'Location', 'SouthEast')
legend boxoff

axes('Position', [UL + 2.8 * boxwidth + 5 * boxspace_h, 1-UL - boxheight, 0.6 * boxwidth, boxheight/2 - boxspace_v/2])

plot([1:50:1000,1000],V([1:50:1000,1000],end-1), 'linewidth', 4, 'marker', 'x', 'markersize', 10, 'color', 'b', 'linestyle', 'none')
axis([0,1000, -0.125,0.125])
set(gca, 'fontsize', 24)
xlabel('Pyramidal cell')
ylabel('Loading Factor')
box on

legend({'2nd P.C.'}, 'Location', 'SouthEast')
legend boxoff

axes('Position', [UL + 3.7 * boxwidth + 2 * boxspace_h, 1-UL - boxheight, 0.6 * boxwidth, boxheight])

hold on

% plot(min(abs(flipud(acos(V' * diff(RE0')' ./ (norm(diff(RE0')))) * 180/pi)), abs(flipud(acos(V' * diff(RE0')' ./ (norm(diff(RE0')))) * 180/pi) - 180)), 'linewidth', 4, 'marker', 'x', 'markersize', 10, 'linestyle', 'none', 'color', 'k')
% hold on
% plot([-10,1010], [90,90], 'color', 'k', 'linewidth', 2, 'linestyle', '--')
% plot([-10,1010], [0,0], 'color', 'k', 'linewidth', 2, 'linestyle', '--')
% 
% set(gca,'fontsize', 24)
% xlabel('Principal Component')
% ylabel('Angle from mean vector (deg.)')
% axis([0.5,10.5,-5,95])

E_Spect_Overlap = flipud(abs(V' * diff(RE0')' ./ (norm(diff(RE0')))));

plot(E_Spect_Overlap(1), 'linewidth', 4, 'marker', 'x', 'markersize', 10, 'linestyle', 'none', 'color', 'r')
plot(2, E_Spect_Overlap(2), 'linewidth', 4, 'marker', 'x', 'markersize', 10, 'linestyle', 'none', 'color', 'b')
plot(3:1000, E_Spect_Overlap(3:end), 'linewidth', 4, 'marker', 'x', 'markersize', 10, 'linestyle', 'none', 'color', 'k')

hold on
plot([-10,1010], [1,1], 'color', 'k', 'linewidth', 2, 'linestyle', '--')
plot([-10,1010], [0,0], 'color', 'k', 'linewidth', 2, 'linestyle', '--')

set(gca,'fontsize', 24)
xlabel('Principal Component')
ylabel('Overlap with Mean Vector')
axis([0.5,10.5,-0.1,1.1])
box on

%% Plot Selectivity Sweep and Decompositino


X = 'SweepEIwidth_squaredSIconvention.mat'; 
load([ROOT, X])

axes('Position', [UL + 0.05 * boxwidth, 1-UL - 2.5 * boxheight, 2.0 * boxwidth, 1.1 * boxheight])

h = imagesc(([1:Nvals] - 26) /20, ([1:Nvals] - 26) /20, sqrt(SItot_E_epsridge') .* (SItot_E_epsridge' > 0))
set(h, 'AlphaData', ~isnan(SItot_E_epsridge'))
set(gca, 'fontsize', 24)
set(gca, 'ydir', 'normal')
set(gca, 'clim', [0,50])
ylabel('Tuning of E to I weights (k_{IE})')
xlabel('Tuning of I to E weights (k_{EI})')
title('Decoder Performance')
axis([(1-26)/20, (Nvals - 26)/20, 0, (Nvals - 26)/20])
h = colorbar;
set(h, 'fontsize', 24)
box on


ROOT     = '/nfs/nhome/live/angus/Documents/Talks and Posters/Cosyne 2018/Further Results/NewParams/Covariance/Tuned I/'

X = 'SweepUntunedtoISO.mat';
load([ROOT, X])

axes('Position', [UL + 1.9 * boxwidth, 1-UL - 2.4 * boxheight, 0.6 * boxwidth, 1.0 * boxheight])

hold on
plot(-0.9:0.1:0, fliplr(NrmSq), 'linewidth', 4, 'color', [0.3,0.3,0.3])
set(gca, 'fontsize', 24)
set(gca, 'xtick', [-1,0,1]);
set(gca, 'xticklabel', {'Iso', 'N.S.', 'Ortho'})
xlabel('Network Type')
ylabel('Signal (||\Delta r||^2)')

X = 'SweepUntunedtoOrtho.mat';
load([ROOT, X])

plot(0:0.1:0.9, NrmSq, 'linewidth', 4, 'color', [0.7,0.7,0.7])
box on
axis([-1,1, -0.1,0.9])

X = 'SweepUntunedtoISO.mat';
load([ROOT, X])

axes('Position', [UL + 2.8 * boxwidth, 1-UL - 2.4 * boxheight, 0.6 * boxwidth, 1.0 * boxheight])

hold on
plot(-0.9:0.1:0, fliplr(TrD), 'linewidth', 4, 'color', [0.3,0.3,0.3])
set(gca, 'fontsize', 24)
set(gca, 'xtick', [-1,0,1]);
set(gca, 'xticklabel', {'Iso', 'N.S.', 'Ortho'})
xlabel('Network Type')
ylabel('Noise (Trace \Sigma)')
box on

X = 'SweepUntunedtoOrtho.mat';
load([ROOT, X])

plot(0:0.1:0.9, TrD, 'linewidth', 4, 'color', [0.7,0.7,0.7])
axis([-1,1, 0,0.02])


axes('Position', [UL + 3.7 * boxwidth, 1-UL - 2.4 * boxheight, 0.6 * boxwidth, 1.0 * boxheight])

X = 'SweepUntunedtoISO.mat';
load([ROOT, X])

hold on
plot(-0.9:0.1:0, fliplr(1./SNRalign), 'linewidth', 4, 'color', [0.3,0.3,0.3])

X = 'SweepUntunedtoOrtho.mat';
load([ROOT, X])

plot(0:0.1:0.9, 1./SNRalign, 'linewidth', 4, 'color', [0.7,0.7,0.7])
set(gca, 'fontsize', 24)
set(gca, 'xtick', [-1,0,1]);
set(gca, 'xticklabel', {'Iso', 'N.S.', 'Ortho'})
xlabel('Network Type')
ylabel('Signal-Noise Alignment')
box on
axis([-1,1, 0,0.025])


%% Plot Sweeps over EI mags

ROOT     = '/nfs/nhome/live/angus/Documents/SSN_Paper/New_Code/'

X1 = 'SweepEImag_nonspecificI_squaredSIconvention.mat'; 
X2 = 'SweepEImag_isoI_squaredSIconvention.mat';
X3 = 'SweepEImag_crossI_squaredSIconvention.mat';

load([ROOT, X1])

axes('Position', [UL + 0.2 * boxwidth, 1-UL - 4.3 * boxheight, 1.2 * boxwidth, 1.2 * boxheight])

h = imagesc(JEI_mean,JIE_mean, log(sqrt(SItot_E .* (SItot_E > 0))))
set(h, 'AlphaData', ~isnan(SItot_E))
set(gca, 'fontsize', 24)
%set(gca, 'clim', [0,4])
set(gca, 'ydir', 'normal')
xlabel('E to I coupling strength (J_{IE})')
ylabel('I to E coupling strength (J_{EI})')
title('Decoder Performance')
h = colorbar;
set(h, 'fontsize', 24)
title(h, 'log SI')

load([ROOT, X2])

axes('Position', [UL + 1.7  * boxwidth, 1-UL - 4.3 * boxheight, 1.2 * boxwidth, 1.2 * boxheight])

h = imagesc(JEI_mean,JIE_mean, log(sqrt(SItot_E .* (SItot_E > 0))))
set(h, 'AlphaData', ~isnan(SItot_E))
set(gca, 'fontsize', 24)
set(gca, 'clim', [0,4])
set(gca, 'ydir', 'normal')
xlabel('E to I coupling strength (J_{IE})')
ylabel('I to E coupling strength (J_{EI})')
title('Decoder Performance')
h = colorbar;
set(h, 'fontsize', 24)
title(h, 'log SI')


load([ROOT, X3])

axes('Position', [UL + 3.2 * boxwidth, 1-UL - 4.3 * boxheight, 1.2 * boxwidth, 1.2 * boxheight])

h = imagesc(JEI_mean,JIE_mean, log(sqrt(SItot_E .* (SItot_E > 0))))
set(h, 'AlphaData', ~isnan(SItot_E))
set(gca, 'fontsize', 24)
set(gca, 'clim', [0,4])
set(gca, 'ydir', 'normal')
xlabel('E to I coupling strength (J_{IE})')
ylabel('I to E coupling strength (J_{EI})')
title('Decoder Performance')
h = colorbar;
set(h, 'fontsize', 24)
title(h, 'log SI')


%% Save Figure 

print('CovarianceFig', '-dsvg')
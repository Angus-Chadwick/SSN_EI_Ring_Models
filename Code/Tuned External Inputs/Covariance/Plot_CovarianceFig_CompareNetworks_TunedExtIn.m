% Plot Activity, covariance and spectra
%% Load Data

clear all

ROOT     = '/mnt/data/angusc/Data/Tuned External Inputs/Covariance/';

X =        'CovarianceData_NSI_TunedExtIn.mat';
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

axis([0.5,10.5, 0,1.2])
set(gca, 'fontsize', 10)
xlabel('Principal Component')
ylabel('Variance Explained (%)')
box on

axes('Position',[UL + 1.9 * boxwidth + 0.6 * boxwidth - 1.25 * inset_size,1-UL - boxheight + boxheight - 1.25 * inset_size, inset_size, inset_size])
plot(flipud(diag(D)) / trace(D) * 100, 'color', 'k', 'linewidth', 4)
set(gca, 'fontsize', 10)
xlabel('Princ. Comp.')
ylabel('Var. Ex. (%)')
axis([-10, 1010, -0.1, 1.2])
box on

axes('Position', [UL + 2.8 * boxwidth, 1-UL - boxheight/2 + boxspace_v/2, 0.6 * boxwidth, boxheight/2 - boxspace_v/2])

plot([1:50:1000,1000],V([1:50:1000,1000],end), 'linewidth', 4, 'marker', 'x', 'markersize', 8, 'color', 'r', 'linestyle', 'none')
axis([0,1000, -0.125,0.125])
set(gca, 'fontsize', 10)
set(gca, 'xticklabel', {})
box on

legend({'1st P.C.'}, 'Location', 'Best')
legend boxoff

axes('Position', [UL + 2.8 * boxwidth, 1-UL - boxheight, 0.6 * boxwidth, boxheight/2 - boxspace_v/2])

plot([1:50:1000,1000],-V([1:50:1000,1000],end-1), 'linewidth', 4, 'marker', 'x', 'markersize', 8, 'color', 'b', 'linestyle', 'none')
axis([0,1000, -0.125,0.125])
set(gca, 'fontsize', 10)
xlabel('Pyramidal cell')
ylabel('Loading Factor')
box on

legend({'2nd P.C.'}, 'Location', 'Best')
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

%% Iso Inhibition

clear all

ROOT     = '/mnt/data/angusc/Data/Tuned External Inputs/Covariance/';

X =        'CovarianceData_IsoI_TunedExtIn_Batch1.mat';
load([ROOT, X])

X =        'CovarianceData_IsoI_TunedExtIn_PooledCovsandmeans.mat'; 
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

axes('Position', [UL, 1-UL - boxheight/2 + boxspace_v/2 - boxheight - 2 * boxspace_v, 0.8 * boxwidth, boxheight/2 - boxspace_v/2])

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

axes('Position', [UL, 1-UL - boxheight  - boxheight - 2 * boxspace_v, 0.8 * boxwidth, boxheight/2 - boxspace_v/2])

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

axes('Position', [UL + 0.8  * boxwidth + 0.5 * boxspace_h,1-UL - boxheight  - boxheight - 2 * boxspace_v, 0.8 * boxwidth, boxheight])

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

axes('Position', [UL + 1.9 * boxwidth, 1-UL - boxheight  - boxheight - 2 * boxspace_v, 0.6 * boxwidth, boxheight])
    
hold on

E_Spect = flipud(diag(D)) / trace(D) * 100;

plot(1,E_Spect(1), 'color', 'r', 'marker', 'x', 'markersize', 8, 'linestyle', 'none', 'linewidth', 4)
plot(2,E_Spect(2), 'color', 'b', 'marker', 'x', 'markersize', 8, 'linestyle', 'none', 'linewidth', 4)
plot(3:1000,E_Spect(3:end), 'color', 'k', 'marker', 'x', 'markersize', 8, 'linestyle', 'none', 'linewidth', 4)

axis([0.5,10.5, 0,1.2])
set(gca, 'fontsize', 10)
xlabel('Principal Component')
ylabel('Variance Explained (%)')
box on

axes('Position',[UL + 1.9 * boxwidth + 0.6 * boxwidth - 1.25 * inset_size,1-UL - boxheight + boxheight - 1.25 * inset_size - boxheight - 2 * boxspace_v, inset_size, inset_size])
plot(flipud(diag(D)) / trace(D) * 100, 'color', 'k', 'linewidth', 4)
set(gca, 'fontsize', 10)
xlabel('Princ. Comp.')
ylabel('Var. Ex. (%)')
axis([-10, 1010, -0.1, 1.2])
box on

axes('Position', [UL + 2.8 * boxwidth, 1-UL - boxheight/2 + boxspace_v/2  - boxheight - 2 * boxspace_v, 0.6 * boxwidth, boxheight/2 - boxspace_v/2])

plot([1:50:1000,1000],-V([1:50:1000,1000],end), 'linewidth', 4, 'marker', 'x', 'markersize', 8, 'color', 'r', 'linestyle', 'none')
axis([0,1000, -0.125,0.125])
set(gca, 'fontsize', 10)
set(gca, 'xticklabel', {})
box on

legend({'1st P.C.'}, 'Location', 'Best')
legend boxoff

axes('Position', [UL + 2.8 * boxwidth, 1-UL - boxheight  - boxheight - 2 * boxspace_v, 0.6 * boxwidth, boxheight/2 - boxspace_v/2])

plot([1:50:1000,1000],-V([1:50:1000,1000],end-1), 'linewidth', 4, 'marker', 'x', 'markersize', 8, 'color', 'b', 'linestyle', 'none')
axis([0,1000, -0.125,0.125])
set(gca, 'fontsize', 10)
xlabel('Pyramidal cell')
ylabel('Loading Factor')
box on

legend({'2nd P.C.'}, 'Location', 'Best')
legend boxoff

axes('Position', [UL + 3.7 * boxwidth, 1-UL - boxheight  - boxheight - 2 * boxspace_v, 0.6 * boxwidth, boxheight])

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

%% Cross inhibition

ROOT     = '/mnt/data/angusc/Data/Tuned External Inputs/Covariance/';

X =        'CovarianceData_CrossI_TunedExtIn.mat';
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

axes('Position', [UL, 1-UL - boxheight/2 + boxspace_v/2 - 2 * boxheight - 4 * boxspace_v, 0.8 * boxwidth, boxheight/2 - boxspace_v/2])

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

axes('Position', [UL, 1-UL - boxheight - 2 * boxheight - 4 * boxspace_v, 0.8 * boxwidth, boxheight/2 - boxspace_v/2])

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

axes('Position', [UL + 0.8  * boxwidth + 0.5 * boxspace_h,1-UL - boxheight  - 2 * boxheight - 4 * boxspace_v, 0.8 * boxwidth, boxheight])

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

axes('Position', [UL + 1.9 * boxwidth, 1-UL - boxheight  - 2 * boxheight - 4 * boxspace_v, 0.6 * boxwidth, boxheight])
    
hold on

E_Spect = flipud(diag(D)) / trace(D) * 100;

plot(1,E_Spect(1), 'color', 'r', 'marker', 'x', 'markersize', 8, 'linestyle', 'none', 'linewidth', 4)
plot(2,E_Spect(2), 'color', 'b', 'marker', 'x', 'markersize', 8, 'linestyle', 'none', 'linewidth', 4)
plot(3:1000,E_Spect(3:end), 'color', 'k', 'marker', 'x', 'markersize', 8, 'linestyle', 'none', 'linewidth', 4)

axis([0.5,10.5, 0,1.2])
set(gca, 'fontsize', 10)
xlabel('Principal Component')
ylabel('Variance Explained (%)')
box on

axes('Position',[UL + 1.9 * boxwidth + 0.6 * boxwidth - 1.25 * inset_size,1-UL - boxheight + boxheight - 1.25 * inset_size - 2 * boxheight - 4 * boxspace_v, inset_size, inset_size])
plot(flipud(diag(D)) / trace(D) * 100, 'color', 'k', 'linewidth', 4)
set(gca, 'fontsize', 10)
xlabel('Princ. Comp.')
ylabel('Var. Ex. (%)')
axis([-10, 1010, -0.1, 1.2])
box on

axes('Position', [UL + 2.8 * boxwidth, 1-UL - boxheight/2 + boxspace_v/2 - 2 * boxheight - 4 * boxspace_v, 0.6 * boxwidth, boxheight/2 - boxspace_v/2])

plot([1:50:1000,1000],V([1:50:1000,1000],end), 'linewidth', 4, 'marker', 'x', 'markersize', 8, 'color', 'r', 'linestyle', 'none')
axis([0,1000, -0.125,0.125])
set(gca, 'fontsize', 10)
set(gca, 'xticklabel', {})
box on

legend({'1st P.C.'}, 'Location', 'Best')
legend boxoff

axes('Position', [UL + 2.8 * boxwidth, 1-UL - boxheight - 2 * boxheight - 4 * boxspace_v, 0.6 * boxwidth, boxheight/2 - boxspace_v/2])

plot([1:50:1000,1000],-V([1:50:1000,1000],end-1), 'linewidth', 4, 'marker', 'x', 'markersize', 8, 'color', 'b', 'linestyle', 'none')
axis([0,1000, -0.125,0.125])
set(gca, 'fontsize', 10)
xlabel('Pyramidal cell')
ylabel('Loading Factor')
box on

legend({'2nd P.C.'}, 'Location', 'Best')
legend boxoff

axes('Position', [UL + 3.7 * boxwidth, 1-UL - boxheight - 2 * boxheight - 4 * boxspace_v, 0.6 * boxwidth, boxheight])

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
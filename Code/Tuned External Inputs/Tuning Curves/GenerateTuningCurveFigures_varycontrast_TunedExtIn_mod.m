% Generate tuning curves, standard devation, and selectivity for three
% networks. Output also the external and recurrent currents received by
% each cell. 

clear all

Match_SI = 0; % choose whether to adjust average E<->I weights to match selectivity across network types 

% simulation parameters

Nloop = 10;  % number of simulations for each parameter set
Nt = 10000;  % number of timesteps

% initialise variables

RE = cell([Nloop, 1]);
RI = RE;
RE_std = RE;
RI_std = RI;
RE_cov = RE;
Rtot_cov = RE;

% fixed network parameters 

JEI_mean = 0.04 ;  
JIE_mean = 0.04 ;

% fixed input parameters

theta_aE = 0;
theta_aI = 0;
noise = 2;
kE_FF = 0.5;
kE_TD = 0;
kI_TD = 0;
IE_TD_area = 0.0;
II_TD_area = 0.0;

stimvals = 2*pi * 180 / 360;

Nstim = 10; 


%%

for m=1:3
    
for q=1:Nstim
    
if m==1 % no tuning

kIE = 0; 
kEI = 0;

if Match_SI == 1

JEI_mean = 0.035;
JIE_mean = 0.035;

end

elseif m==2 % iso-orientation inhibition

kIE = 0.5; 
kEI = 0.5;

%kEI = 0.3;
%kIE = 0.3;

if Match_SI == 1

JEI_mean = 0.0265;
JIE_mean = 0.0265;

%JEI_mean = 0.0315;
%JIE_mean = 0.0315;

end

elseif m == 3 % cross-orientation inhibition

kIE = 0.5;
kEI = -0.5;

%kIE = 0.3;
%kEI = -0.3;

if Match_SI == 1

JEI_mean = 0.04;
JIE_mean = 0.04;

%JEI_mean = 0.0375;
%JIE_mean = 0.0375;

end

end
    
network = create_network(kEI,kIE,JEI_mean,JIE_mean);

NE = network.cells.NE;
NI = network.cells.NI;

% create inputs

theta_s = stimvals;
IE_FF_area(q) = 0.05 * q;

inputs  = create_inputs(theta_s, theta_aE, theta_aI, noise, kE_FF, IE_FF_area(q), kE_TD, IE_TD_area, kI_TD, II_TD_area, network);

 
%% simulate

NoiseModel = 'Add';


parfor n=1:Nloop

    [rE, rI]       = SimulateNetwork_mod(network, inputs, Nt, NoiseModel);
    
    RE{n}(:,q) = mean(rE(:,300:end),2);
    RI{n}(:,q) = mean(rI(:,300:end),2);
    
    RE_std{n}(:,q) = std(rE(:,300:end),[],2);
    RI_std{n}(:,q) = std(rI(:,300:end),[],2);

    RE_cov{n}(:,:,q) = cov([rE(:,300:end); rI(:,300:end)]')

end

    



RE0 = mean(cat(3,RE{:}),3);
RI0 = mean(cat(3,RI{:}),3);

RE0_std = mean(cat(3,RE_std{:}),3);
RI0_std = mean(cat(3,RI_std{:}),3);

RE0_cov = mean(cat(4,RE_cov{:}),4);

        RE0_std_full{m,q} = squeeze(RE0_std(:,q));        
        RE0_full{m,q} = squeeze(RE0(:,q));
        
        RI0_std_full{m,q} = squeeze(RI0_std(:,q));
        RI0_full{m,q} = squeeze(RI0(:,q));
             
end
end

MeanRateE = cellfun(@nanmean, RE0_full);
MeanRateI = cellfun(@nanmean, RI0_full);

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

theta_pE = inputs.theta_pE;
theta_pI = inputs.theta_pI;

figure
subplot(2,3,1)
hold on
set(gca, 'fontsize', 24)
ylabel('Mean Response')
title('Pyramidal Cells')

for q=1:Nstim
    hold on
    plot(theta_pE * 180/pi, RE0_full{1,q}, 'color', Xr, 'linewidth', 1, 'linestyle', '-')
end

set(gca, 'xtick', 0:90:360)
set(gca, 'xticklabel', {})
box on
subplot(2,3,4)
hold on
set(gca, 'fontsize', 24)
ylabel('Mean Response')
title('Interneurons')
for q=1:Nstim
    hold on
    plot(theta_pI * 180/pi, RI0_full{1,q}, 'color', Xr, 'linewidth', 1, 'linestyle', '-')
end
set(gca, 'xtick', 0:90:360)
%set(gca, 'xticklabel', {})
box on
subplot(2,3,2)
hold on
set(gca, 'fontsize', 24)
for q=1:Nstim
    hold on
    plot(theta_pE * 180/pi, RE0_full{2,q}, 'color', Xr, 'linewidth', 1, 'linestyle', '-')
end
set(gca, 'xtick', 0:90:360)
set(gca, 'xticklabel', {})
box on
subplot(2,3,5)
hold on
set(gca, 'fontsize', 24)
for q=1:Nstim
    hold on
    plot(theta_pI * 180/pi, RI0_full{2,q}, 'color', Xr, 'linewidth', 1, 'linestyle', '-')
end
set(gca, 'xtick', 0:90:360)
%set(gca, 'xticklabel', {})
box on
subplot(2,3,3)
hold on
set(gca, 'fontsize', 24)
for q=1:Nstim
    hold on
    plot(theta_pE * 180/pi, RE0_full{3,q}, 'color', Xr, 'linewidth', 1, 'linestyle', '-')
end
subplot(2,3,6)
hold on
set(gca, 'fontsize', 24)
ylabel('Ratio of Responses')
for q=1:Nstim
    hold on
    plot(theta_pI * 180/pi, RI0_full{3,q}, 'color', Xr, 'linewidth', 1, 'linestyle', '-')
end

figure
subplot(2,1,1)
plot(IE_FF_area, MeanRateE')
subplot(2,1,2)
plot(IE_FF_area, MeanRateI')
xlabel('Stimulus Contrast')
ylabel('Mean Response')
legend('Nonspecific I', 'Iso I', 'Cross I')


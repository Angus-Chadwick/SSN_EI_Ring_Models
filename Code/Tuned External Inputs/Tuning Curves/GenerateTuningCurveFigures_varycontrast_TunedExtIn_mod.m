% Generate tuning curves, standard devation, and selectivity for three
% networks. Output also the external and recurrent currents received by
% each cell. 

clear all

Match_SI = 1 % choose whether to adjust average E<->I weights to match selectivity across network types 

% simulation parameters

Nloop = 30;  % number of simulations for each parameter set
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

Nstim = 5; 


%%

for m=3:3
    
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
kEI= 0.5;

if Match_SI == 1

 JEI_mean = 0.0265;
 JIE_mean = 0.0265;

end

elseif m == 3 % cross-orientation inhibition

kIE = 0.5;
kEI = -0.5;


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
IE_FF_area(q) = 10 * q;

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

theta_pE = inputs.theta_pE;
theta_pI = inputs.theta_pI;

figure

subplot(3,2,1)
plot(IE_FF_area, MeanRateE', 'linewidth', 4)
set(gca, 'fontsize', 16)
xlabel('Stimulus Contrast')
ylabel('Mean Response')
legend('Nonspecific I', 'Iso I', 'Cross I')
title('Pyramidal')
subplot(3,2,2)
plot(IE_FF_area, MeanRateI', 'linewidth', 4)
xlabel('Stimulus Contrast')
ylabel('Mean Response')
legend('Nonspecific I', 'Iso I', 'Cross I')
title('Interneuron')
set(gca, 'fontsize', 16)


cmap = jet(Nstim);

for m=1:3
for q=20:10:70
% pyr TC
subplot(3,3,m+3)
hold on
set(gca, 'fontsize', 14)
ylabel('Mean Response')
    
    hold on
    plot(theta_pE * 180/pi, RE0_full{m,q}, 'color', cmap(q,:), 'linewidth', 3, 'linestyle', '-')

axis([0,360, 0,0.08])
box on
    
% pyr TC ratio
subplot(3,3,m+6)
hold on
set(gca, 'fontsize', 14)
ylabel('Response Gain')
box on
    hold on
    plot(theta_pE * 180/pi, RE0_full{m,q} ./ RE0_full{m,50}, 'color', cmap(q,:), 'linewidth', 3, 'linestyle', '-')

set(gca, 'xtick', 0:90:360)
%set(gca, 'xticklabel', {})
box on
axis([0,360, 0,3])
% 
% % int TC
% subplot(4,3,m+6)
% hold on
% set(gca, 'fontsize', 14)
% ylabel('Mean Response')
% title('Interneurons')
% 
%     hold on
%     plot(theta_pI * 180/pi, RI0_full{m,q}, 'color', cmap(q,:), 'linewidth', 1, 'linestyle', '-')
% 
% set(gca, 'xtick', 0:90:360)
% box on
% axis([0,360, 0,0.25])
% 
% % int TC ratio 
% subplot(4,3,m+9)
% hold on
% set(gca, 'fontsize', 14)
% ylabel('Response Gain')
% title('Interneurons')
% 
%     hold on
%     plot(theta_pI * 180/pi, RI0_full{m,q} ./ RI0_full{m,Nstim/2}, 'color', cmap(q,:), 'linewidth', 1, 'linestyle', '-')
%     
% set(gca, 'xtick', 0:90:360)
% box on
% axis([0,360, 0,4])

end
end



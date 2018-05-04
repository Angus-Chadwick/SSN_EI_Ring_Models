%% Ring model with attention

clear all

% simulation parameters

Nloop = 10;  % number of simulations for each parameter set
Nvals = 50;  % number of parameter sets
Nt = 10000;  % number of timesteps

% fixed network parameters 

kIE = 0.4; 
kEI = 0.4;
JEI_mean = 0.026;  
JIE_mean = 0.026;

network = create_network(kEI,kIE,JEI_mean,JIE_mean);

% fixed input parameters

stimvals = 2*pi * [160,200] / 360;
theta_aE = (stimvals(1) + stimvals(2))/2 + pi;
theta_aI = (stimvals(1) + stimvals(2))/2;
noise = 2;
kE_FF = 0.5;
kE_TD = 1;
kI_TD = 1;
IE_FF_area = 0.005 * 100;
Nstim = length(stimvals);

% simulate with different top-down inputs

for p1 = 1:Nvals
for p2 = 1:Nvals
for q=1:Nstim

theta_s = stimvals(q);

IE_TD_area = 0.02 * p1;
II_TD_area = 0.02 * p2;

inputs  = create_inputs(theta_s, theta_aE, theta_aI, noise, kE_FF, IE_FF_area, kE_TD, IE_TD_area, kI_TD, II_TD_area, network);

%% simulate

NoiseModel = 'Add';

parfor n=1:Nloop

    [rE_TD, rI_TD] = SimulateNetwork_mod(network, inputs, Nt, NoiseModel);

    RE_TD{n}(:,q) = mean(rE_TD(:,300:end),2);
    RI_TD{n}(:,q) = mean(rI_TD(:,300:end),2);
      
    RE_std_TD{n}(:,q) = std(rE_TD(:,300:end),[],2);
    RI_std_TD{n}(:,q) = std(rI_TD(:,300:end),[],2);
        
    RE_cov_TD{n}(:,:,q) = cov(rE_TD(:,300:end)');
    Rtot_cov_TD{n}(:,:,q) = cov([rE_TD(:,300:end); rI_TD(:,300:end)]')
    
end

end

RE0_TD = mean(cat(3,RE_TD{:}),3);
RI0_TD = mean(cat(3,RI_TD{:}),3);

RE0_std_TD = mean(cat(3,RE_std_TD{:}),3);
RI0_std_TD = mean(cat(3,RI_std_TD{:}),3);

RE0_cov_TD = mean(cat(4,RE_cov_TD{:}),4);
Rtot0_cov_TD = mean(cat(4,Rtot_cov_TD{:}),4);

RE_covtot1_TD{p1,p2} = RE0_cov_TD(:,:,1);
RE_covtot2_TD{p1,p2} = RE0_cov_TD(:,:,2);
SItot_E_TD(p1,p2) = squeeze(RE0_TD(:,1) - RE0_TD(:,2))' * inv(0.5 * (RE_covtot1_TD{p1,p2}+ RE_covtot2_TD{p1,p2})) * squeeze(RE0_TD(:,1) - RE0_TD(:,2)); 
SItot_E_ind_TD(p1,p2) =  squeeze(RE0_TD(:,1) - RE0_TD(:,2))' * inv(0.5 * diag(diag(RE_covtot1_TD{p1,p2} + RE_covtot2_TD{p1,p2}))) * squeeze(RE0_TD(:,1) - RE0_TD(:,2)); 

end
end
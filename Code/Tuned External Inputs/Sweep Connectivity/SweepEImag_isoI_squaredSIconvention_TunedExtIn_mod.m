%% Ring model with attention

clear all

% simulation parameters

Nloop = 5;  % number of simulations for each parameter set
Nvals = 40;  % number of parameter sets
Nt = 10000;  % number of timesteps

% initialise variables

RE_covtot = cell(Nvals);
RE = cell([Nloop, 1]);
RI = RE;
RE_std = RE;
RI_std = RI;
RE_cov = RE;
Rtot_cov = RE;

% fixed input parameters

stimvals = 2*pi * [160,200] / 360;
Nstim = length(stimvals);
theta_aE = 0;
theta_aI = 0;
noise = 2;
kE_FF = 0.5;
kE_TD = 0;
kI_TD = 0;
IE_FF_area = 0.005 * 100;
IE_TD_area = 0.0;
II_TD_area = 0.0;

% fixed connectivity parameters

kEI = 0.5;
kIE = 0.5;

for nI = 1:Nvals
    nI
    for nE = 1:Nvals
    
    for q=1:Nstim
 
% create network        


JEI_mean(nI,nE) = 0.1/Nvals * nI;
JIE_mean(nI,nE) = 0.1/Nvals * nE;

network = create_network(kEI,kIE,JEI_mean(nI,nE),JIE_mean(nI,nE));

NE = network.cells.NE;
NI = network.cells.NI;

% create inputs

theta_s = stimvals(q);

inputs  = create_inputs(theta_s, theta_aE, theta_aI, noise, kE_FF, IE_FF_area, kE_TD, IE_TD_area, kI_TD, II_TD_area, network);

    

%% simulate

NoiseModel = 'Add';


parfor n=1:Nloop

            [rE, rI]       = SimulateNetwork_mod(network, inputs, Nt, NoiseModel);

    
    RE{n}(:,q) = mean(rE(:,300:end),2);
    RI{n}(:,q) = mean(rI(:,300:end),2);
    
    RE_std{n}(:,q) = std(rE(:,300:end),[],2);
    RI_std{n}(:,q) = std(rI(:,300:end),[],2);

    RE_cov{n}(:,:,q) = cov(rE(:,300:end)');
    Rtot_cov{n}(:,:,q) = cov([rE(:,300:end); rI(:,300:end)]')
    

end

    end

%% analyse outputs


RE0 = mean(cat(3,RE{:}),3);
RI0 = mean(cat(3,RI{:}),3);

RE0_std = mean(cat(3,RE_std{:}),3);
RI0_std = mean(cat(3,RI_std{:}),3);

RE0_cov = mean(cat(4,RE_cov{:}),4);
Rtot0_cov = mean(cat(4,Rtot_cov{:}),4);


mmE = zeros(Nstim);
mmI = zeros(Nstim);
m0E = zeros(Nstim);
m0I = zeros(Nstim);

iv1 = 1;iv2 = 2;

        RE0_std_vert = squeeze(RE0_std(:,iv1));
        RE0_std_ang = squeeze(RE0_std(:,iv2));

        RE0_vert = squeeze(RE0(:,iv1));
        RE0_ang = squeeze(RE0(:,iv2));
     
        RI0_std_vert = squeeze(RI0_std(:,iv1));
        RI0_std_ang = squeeze(RI0_std(:,iv2));

        RI0_vert = squeeze(RI0(:,iv1));
        RI0_ang  = squeeze(RI0(:,iv2));

        % E selectivity
        
        dmu_0E = RE0_ang - RE0_vert;
        
        poolvar_0E = (0.5* (RE0_std_ang.^2 + RE0_std_vert.^2));

        SI0E = dmu_0E.^2 ./ poolvar_0E;

        m0E(iv1,iv2) = nansum(SI0E);

        % I selectivity
        
        dmu_0I = RI0_ang - RI0_vert;
        
        poolvar_0I = (0.5* (RI0_std_ang.^2 + RI0_std_vert.^2));

        SI0I = dmu_0I.^2 ./ poolvar_0I;

        m0I(iv1,iv2) = nansum(SI0I);

RE_covtot1{nI,nE} = RE0_cov(:,:,1);
RE_covtot2{nI,nE} = RE0_cov(:,:,2);
SItot_E(nI,nE) = squeeze(RE0(:,1) - RE0(:,2))' * inv(0.5 * (RE_covtot1{nI,nE} + RE_covtot2{nI,nE})) * squeeze(RE0(:,1) - RE0(:,2)); 
SItot_E_ind(nI,nE) =  squeeze(RE0(:,1) - RE0(:,2))' * inv(0.5 * diag(diag(RE_covtot1{nI,nE} + RE_covtot2{nI,nE}))) * squeeze(RE0(:,1) - RE0(:,2)); 
SI_I(nI,nE) = m0I(1,2) / NI;
SI_E(nI,nE) = m0E(1,2) / NE;



%% Condition Covariance Matrix

if ~isnan(sum(RE_covtot1{nI,nE} + RE_covtot2{nI,nE}))

    SItot_E_pseudo(nI,nE) = squeeze(RE0(:,1) - RE0(:,2))' * pinv(0.5 * (RE_covtot1{nI,nE} + RE_covtot2{nI,nE})) * squeeze(RE0(:,1) - RE0(:,2)); 

else
    
    SItot_E_pseudo(nI,nE) = nan;
    
end

SItot_E_epsridge(nI,nE) = squeeze(RE0(:,1) - RE0(:,2))' * inv(0.5 * (RE_covtot1{nI,nE} + RE_covtot2{nI,nE} + 2 * eps * eye(1000))) * squeeze(RE0(:,1) - RE0(:,2)); 

    
    end
end


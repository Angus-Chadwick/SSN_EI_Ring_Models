
clear all

% simulation parameters

Nt = 10000;  % number of timesteps
Nloop = 100;
Ncon = 11;

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
IE_FF_area = 0.005 * 100;
IE_TD_area = 0.0;
II_TD_area = 0.0;


stimvals = 2*pi * [160, 200] / 360;
Nstim = length(stimvals); 


for N = 1:Ncon

kIE = 0.1 * (N-1); 
kEI = 0.1 * (N-1);
    
network = create_network(kEI,kIE,JEI_mean,JIE_mean);

NE = network.cells.NE;
NI = network.cells.NI;

for q=1:Nstim

%% create inputs

theta_s = stimvals(q);

inputs  = create_inputs(theta_s, theta_aE, theta_aI, noise, kE_FF, IE_FF_area, kE_TD, IE_TD_area, kI_TD, II_TD_area, network);


%% simulate

NoiseModel = 'Add';

parfor n=1:Nloop

    [rE, rI] = SimulateNetwork_mod(network, inputs, Nt, NoiseModel);
    
    RE{n}(:,q) = mean(rE(:,300:end),2);
    RI{n}(:,q) = mean(rI(:,300:end),2);

    RE_std{n}(:,q) = std(rE(:,300:end),[],2);
    RI_std{n}(:,q) = std(rI(:,300:end),[],2);

    RE_cov{n}(:,:,q) = cov(rE(:,300:end)');
    Rtot_cov{n}(:,:,q) = cov([rE(:,300:end); rI(:,300:end)]')
    
end

end

%% Analyse response mean, covariance and selectivity

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

for iv1 = 1:Nstim
    for iv2 = 1:Nstim

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
        

    end
end

RE_covtot1 = RE0_cov(:,:,1);
RE_covtot2 = RE0_cov(:,:,2);
SItot_E = squeeze(RE0(:,1) - RE0(:,2))' * inv(0.5 * (RE_covtot1 + RE_covtot2)) * squeeze(RE0(:,1) - RE0(:,2)); 
SItot_E_ind =  squeeze(RE0(:,1) - RE0(:,2))' * inv(0.5 * diag(diag(RE_covtot1 + RE_covtot2))) * squeeze(RE0(:,1) - RE0(:,2)); 
SI_I = m0I(1,2);
SI_E = m0E(1,2);

%% Sample mean and covariance

[V,D] = (eig(0.5 * (RE_covtot1 + RE_covtot2 )));

for i=1:size(D,1)
    cos_angle(i) = (RE0(:,1) - RE0(:,2))' * V(:,i) / norm((RE0(:,1) - RE0(:,2))');  % projection of selectivity into eigenmode
end



TrD(N) = trace(D);
NrmSq(N) = norm(RE0(:,1) - RE0(:,2))^2;
SNRalign(N) = sum(cos_angle.^2 ./ ( diag(D) / trace(D))'); 
SIpop(N) = SItot_E;
SIind(N) = SItot_E_ind;

%% Ridged covariance

[V,D] = (eig(0.5 * (RE_covtot1 + RE_covtot2 + 2 * eps * diag(diag(length(RE_covtot1))))));

for i=1:size(D,1)
    cos_angle_epsridge(i) = (RE0(:,1) - RE0(:,2))' * V(:,i) / norm((RE0(:,1) - RE0(:,2))');  % projection of selectivity into eigenmode
end



TrD_epsridge(N) = trace(D);
NrmSq_epsridge(N) = norm(RE0(:,1) - RE0(:,2))^2;
SNRalign_epsridge(N) = sum(cos_angle_epsridge.^2 ./ ( diag(D) / trace(D))'); 

N % print progress

end

% Generate tuning curves, standard devation, and selectivity for three
% networks. Output also the external and recurrent currents received by
% each cell. 

clear all

Match_SI = 1; % choose whether to adjust average E<->I weights to match selectivity across network types 

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
IE_FF_area = 0.005 * 100;
IE_TD_area = 0.0;
II_TD_area = 0.0;

stimvals = 2*pi * [160, 200] / 360;
Nstim = length(stimvals); 
       
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


    RE_cov{n}(:,:,q) = cov([rE(:,300:end); rI(:,300:end)]')

end

   


    end

% Analyse results

RE0 = mean(cat(3,RE{:}),3);
RI0 = mean(cat(3,RI{:}),3);

RE0_std = mean(cat(3,RE_std{:}),3);
RI0_std = mean(cat(3,RI_std{:}),3);

RE0_cov = mean(cat(4,RE_cov{:}),4);


mmE = zeros(Nstim);
mmI = zeros(Nstim);
m0E = zeros(Nstim);
m0I = zeros(Nstim);

iv1 = 1; iv2 = 2;

        RE0_std_vert{m} = squeeze(RE0_std(:,iv1));
        RE0_std_ang{m} = squeeze(RE0_std(:,iv2));

        RE0_vert{m} = squeeze(RE0(:,iv1));
        RE0_ang{m} = squeeze(RE0(:,iv2));
     
        RI0_std_vert{m} = squeeze(RI0_std(:,iv1));
        RI0_std_ang{m} = squeeze(RI0_std(:,iv2));

        RI0_vert{m} = squeeze(RI0(:,iv1));
        RI0_ang{m} = squeeze(RI0(:,iv2));

        % E selectivity
        
        dmu_0E{m} = RE0_ang{m} - RE0_vert{m};
        
        poolvar_0E{m} = (0.5* (RE0_std_ang{m}.^2 + RE0_std_vert{m}.^2));

        SI0E{m} = dmu_0E{m} ./ sqrt(poolvar_0E{m});
        
        % I selectivity
        
        dmu_0I{m} = RI0_ang{m} - RI0_vert{m};
        
        poolvar_0I{m} = (0.5 * (RI0_std_ang{m}.^2 + RI0_std_vert{m}.^2));

        SI0I{m} = dmu_0I{m} ./ sqrt(poolvar_0I{m});

        RE_covtot1{m} = RE0_cov(:,:,1);
        RE_covtot2{m} = RE0_cov(:,:,2);
       
        ExternalE_Input{m} = inputs.IE_FF;
        RecurrentEE_Input{m} = network.connectivity.JEE * RE0_ang{m};
        RecurrentEI_Input{m} = network.connectivity.JEI * RI0_ang{m};
    
        ExternalI_Input{m} = inputs.II_FF;
        RecurrentIE_Input{m} = network.connectivity.JIE * RE0_ang{m};
        RecurrentII_Input{m} = network.connectivity.JII * RI0_ang{m};
        
end





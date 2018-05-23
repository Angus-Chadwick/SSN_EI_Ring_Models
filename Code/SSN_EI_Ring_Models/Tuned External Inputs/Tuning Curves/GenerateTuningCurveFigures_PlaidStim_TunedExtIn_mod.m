
% Generate tuning curves when multiple stimuli are presented simultaneously

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
kE_FF = 2;
kE_TD = 0;
kI_TD = 0;
IE_FF_area = 0.005 * 100;
IE_TD_area = 0.0;
II_TD_area = 0.0;

stimvals = 2*pi * [90, 270] / 360;
Nstim = length(stimvals); 
       
%%
%%

for m=1:3

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

end

end


network = create_network(kEI,kIE,JEI_mean,JIE_mean);

NE = network.cells.NE;
NI = network.cells.NI;

inputs1  = create_inputs(stimvals(1), theta_aE, theta_aI, noise, kE_FF, IE_FF_area, kE_TD, IE_TD_area, kI_TD, II_TD_area, network);
inputs2  = create_inputs(stimvals(2), theta_aE, theta_aI, noise, kE_FF, IE_FF_area, kE_TD, IE_TD_area, kI_TD, II_TD_area, network);
inputs  = inputs1;
inputs.IE_FF = (inputs1.IE_FF + inputs2.IE_FF);

%% simulate

NoiseModel = 'Add';

parfor n=1:Nloop

            [rE, rI]       = SimulateNetwork_mod(network, inputs, Nt, NoiseModel);

    
    RE{m,n} = mean(rE(:,3000:end),2);
    RI{m,n} = mean(rI(:,3000:end),2);
    
    RE_std{m,n} = std(rE(:,3000:end),[],2);
    RI_std{m,n} = std(rI(:,3000:end),[],2);


    RE_cov{m,n} = cov([rE(:,3000:end); rI(:,3000:end)]')

end

% Analyse results

RE0{m} = mean(cat(3,RE{:}),3);
RI0{m} = mean(cat(3,RI{:}),3);

RE0_std{m} = mean(cat(3,RE_std{:}),3);
RI0_std{m} = mean(cat(3,RI_std{:}),3);

RE0_cov{m} = mean(cat(4,RE_cov{:}),4);

       
        ExternalE_Input{m} = inputs.IE_FF;
        RecurrentEE_Input{m} = network.connectivity.JEE * RE0{m};
        RecurrentEI_Input{m} = network.connectivity.JEI * RI0{m};
    
        ExternalI_Input{m} = inputs.II_FF;
        RecurrentIE_Input{m} = network.connectivity.JIE * RE0{m};
        RecurrentII_Input{m} = network.connectivity.JII * RI0{m};
        

end






clear all

% simulation parameters

Nloop = 10;  % number of simulations for each parameter set
Nvals = 51;
Nt = 5000;  % number of timesteps

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
noise = 0;
kE_FF = 0.5;
kE_TD = 0;
kI_TD = 0;
IE_FF_area = 0.005 * 100;
IE_TD_area = 0.0;
II_TD_area = 0.0;
JEI_mean = 0.04;
JIE_mean = 0.04;

kEI = 0.5;

stepEE = 25;
stepIE = 50;

for nEI = 1:Nvals
    nEI
    for nIE = 1:(5*Nvals)
        
q=1;

% create network        

kEE = (nIE-1) / stepEE; % E to E concentration
kIE = (nEI-1) /stepIE; % I to E concentration

NE = 1000;
network = create_network_varyEE(kEI,kIE,JEI_mean,JIE_mean, 15/NE * besseli(0,1), kEE);

NE = network.cells.NE;
NI = network.cells.NI;

% create inputs

theta_s = stimvals(q);

inputs  = create_inputs(theta_s, theta_aE, theta_aI, noise, kE_FF, IE_FF_area, kE_TD, IE_TD_area, kI_TD, II_TD_area, network);

   
% %% simulate
% 
 NoiseModel = 'Add'; 

[rE, rI]       = SimulateNetwork_mod(network, inputs, Nt, NoiseModel);

R0 = [mean(rE(:,(Nt/2):end),2)', mean(rI(:,(Nt/2):end),2)'];

if ~isnan(sum(R0))
        W = [network.connectivity.JEE, -network.connectivity.JEI; network.connectivity.JIE, -network.connectivity.JII];
        Phip = diag(2 * R0.^(1/2));
        T = diag([network.cells.tauE * ones(NE,1); network.cells.tauI * ones(NI,1)]);
        Jstar = inv(T) * (Phip * W - eye(NE+NI));
        
        [V,D] = eig(Jstar);
        Evals{nEI,nIE} = diag(D);
end

    end
end


kEE = ([1:nIE]-1) / stepEE;
kIE = ([1:nIE]-1) / stepIE;

for i=1:size(Evals,1)
for j=1:size(Evals,2)
if numel(Evals{i,j}) > 0
l(i,j) = max(real(Evals{i,j}));
else l(i,j) = nan;
end
end
end


imagesc(kEE, kIE, (l < 0) + 2 * isnan(l))
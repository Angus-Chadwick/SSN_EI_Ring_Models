%% Learn Optimal representation of stimuli

clear all

% fixed input parameters
NE = 1000;

stimvals = 2*pi * [160,200] / 360;
Nstim = length(stimvals);

noise = 0;
kE_FF = 2; % standard value = 0.5
kI_FF = kE_FF;
IE_FF_area = 0.005 * 100;
II_FF_area = IE_FF_area / 5;
JEE_mean = 15/NE * besseli(0,1);
JEI_mean = 0.04;  % standard values
JIE_mean = 0.04;
JII_mean = JEE_mean * 1.1;
kEE = 0;
kEI = 0;
kIE = 0;
kII = 0;
JEE_mean = 0;
JEI_mean = 0;
JIE_mean = 0;
JII_mean = 0;

network = create_network_varyall(kEE,kEI,kIE,kII, JEE_mean, JEI_mean, JIE_mean, JII_mean);

NE = network.cells.NE;
NI = network.cells.NI;

%% create inputs

theta_s = stimvals(1);

inputs  = create_inputs_varyall(theta_s, noise, kE_FF, kI_FF, IE_FF_area, II_FF_area, network);
inputvec = [inputs.IE_FF', inputs.II_FF'];
inputvecout = inputvec.^2;
%% Learn Weights

W = zeros(NE+NI);
learningrate = 0.001;

signs =  [ones(NE), -ones([NE,NI]); ones([NI,NE]), -ones(NI)]; 

for i=1:1000
    
    gradient = (W .* (inputvecout' * inputvecout) - 1/4 * inputvec' * inputvecout);
    W = W - learningrate * gradient;
    W(W * signs < 0) = 0;
end
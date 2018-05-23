%% Create Inputs 

% input parameters: stimvals, theta_aE, theta_aI, noise 

function inputs = create_inputs(theta_s, theta_aE, theta_aI, noise, kE_FF, IE_FF_area, kE_TD, IE_TD_area, kI_TD, II_TD_area, network)

NE = network.cells.NE;
NI = network.cells.NI; 

% External Inputs

%stimvals = 2*pi * [160,200] / 360;
       
theta_pE = linspace(0, 2*pi, NE+1);
theta_pE = theta_pE(1:(end-1)); % for circular invariance
theta_pI = linspace(0, 2*pi, NI+1);
theta_pI = theta_pI(1:(end-1));

[~, id1]  = min(abs(theta_s - theta_pI));  % centre stimuli on a cell (centred both on an E and an I due to integer NE/NI)
theta_s = theta_pI(id1);

[~, id1]  = min(abs(theta_aE - theta_pI));  % centre stimuli on a cell (centred both on an E and an I due to integer NE/NI)
[~, id2]  = min(abs(theta_aI - theta_pI));
theta_aE = theta_pI(id1);
theta_aI = theta_pI(id2);

% Stimulus drive

IE_FF = (IE_FF_area / (2*pi* besseli(0,kE_FF))) * exp(kE_FF * cos(theta_pE - theta_s))';  
II_FF = -0 * ones([NI,1]); 
    
IE_TD = IE_TD_area / (2*pi* besseli(0,kE_TD)) * exp(kE_TD * cos(theta_pE - theta_aE))';
II_TD = II_TD_area / (2*pi* besseli(0,kE_TD)) * exp(kI_TD * cos(theta_pI - theta_aI))';
    

% create output struct

inputs = struct;
inputs.IE_FF = IE_FF;
inputs.II_FF = II_FF;
inputs.IE_TD = IE_TD;
inputs.II_TD = II_TD;
inputs.noise = noise;
inputs.theta_pE = theta_pE;
inputs.theta_pI = theta_pI;

end
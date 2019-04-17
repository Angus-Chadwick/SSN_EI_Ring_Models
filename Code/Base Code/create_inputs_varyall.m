%% Create Inputs to E and I

function inputs = create_inputs_varyall(theta_s, noise, kE_FF, kI_FF, IE_FF_area, II_FF_area, network)

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


% Stimulus drive

IE_FF = (IE_FF_area / (2*pi* besseli(0,kE_FF))) * exp(kE_FF * cos(theta_pE - theta_s))';  
II_FF = (II_FF_area / (2*pi* besseli(0,kI_FF))) * exp(kI_FF * cos(theta_pI - theta_s))';  
    
   
% create output struct

inputs = struct;
inputs.IE_FF = IE_FF;
inputs.II_FF = II_FF;
inputs.noise = noise;
inputs.theta_pE = theta_pE;
inputs.theta_pI = theta_pI;

end
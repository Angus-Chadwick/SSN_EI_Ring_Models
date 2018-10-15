function [rE, rI] = SimulateNetwork_mod(network, inputs, Nt, NoiseModel)

NE = network.cells.NE;
NI = network.cells.NI;

rE = zeros([NE, Nt]);
rI = zeros([NI, Nt]);

IE_FF = inputs.IE_FF;
II_FF = inputs.II_FF;
IE_TD = inputs.IE_TD;
II_TD = inputs.II_TD;

noise = inputs.noise;

JEE = network.connectivity.JEE;
JEI = network.connectivity.JEI;
JIE = network.connectivity.JIE;
JII = network.connectivity.JII;

tauE = network.cells.tauE;
tauI = network.cells.tauI;
gamma = network.cells.gamma;

if strcmp(NoiseModel, 'Mult')

    for t=1:Nt

        nu_E_FF = noise * randn(NE,1);
        nu_E_FB = noise * randn(NE,1);
        nu_I = noise * randn(NI,1);

        dumE = JEE * rE(:,t) - JEI * rI(:,t) + IE_FF.*(1+nu_E_FF) + IE_TD.*(1+nu_E_FB);
        dumI = JIE * rE(:,t) - JII * rI(:,t) + II_FF.*(1+nu_I) + II_TD.*(1+nu_I);

        drE = (-rE(:,t) + (dumE .* (dumE > 0 )).^gamma) / tauE;
        drI = (-rI(:,t) + (dumI .* (dumI > 0)).^gamma) / tauI;

        rE(:,t+1) = rE(:,t) + drE;
        rI(:,t+1) = rI(:,t) + drI;

    end

elseif strcmp(NoiseModel, 'Add')
    
    noiseE = noise * mean(IE_FF);
    noiseI = noiseE / 2;

    for t=1:Nt

        nu_E_FF = noiseE * randn(NE,1);
        nu_I_FF = noiseI * randn(NI,1);

        dumE = JEE * rE(:,t) - JEI * rI(:,t) + IE_FF + IE_TD + nu_E_FF;
        dumI = JIE * rE(:,t) - JII * rI(:,t) + II_FF + II_TD + nu_I_FF;

        drE = (-rE(:,t) + (dumE .* (dumE > 0 )).^gamma) / tauE;
        drI = (-rI(:,t) + (dumI .* (dumI > 0)).^gamma) / tauI;

        rE(:,t+1) = rE(:,t) + drE;
        rI(:,t+1) = rI(:,t) + drI;

    end

end
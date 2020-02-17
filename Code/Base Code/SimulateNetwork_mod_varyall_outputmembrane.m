function [rE, rI,membraneE,membraneI] = SimulateNetwork_mod_varyall_outputmembrane(network, inputs, Nt, NoiseModel,R0)

NE = network.cells.NE;
NI = network.cells.NI;

rE = zeros([NE, Nt]);
rE(:,1) = R0(1:1000);

rI = zeros([NI, Nt]);
rI(:,1) = R0(1001:end);

membraneE = zeros([NE, Nt]);
membraneI = zeros([NI, Nt]);

IE_FF = inputs.IE_FF;
II_FF = inputs.II_FF;

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
        nu_I = noise * randn(NI,1);

        dumE = JEE * rE(:,t) - JEI * rI(:,t) + IE_FF.*(1+nu_E_FF) ;
        dumI = JIE * rE(:,t) - JII * rI(:,t) + II_FF.*(1+nu_I);

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

        dumE = JEE * rE(:,t) - JEI * rI(:,t) + IE_FF  + nu_E_FF;
        dumI = JIE * rE(:,t) - JII * rI(:,t) + II_FF  + nu_I_FF;

        drE = (-rE(:,t) + (dumE .* (dumE > 0 )).^gamma) / tauE;
        drI = (-rI(:,t) + (dumI .* (dumI > 0)).^gamma) / tauI;

        rE(:,t+1) = rE(:,t) + drE;
        rI(:,t+1) = rI(:,t) + drI;
        membraneE(:,t) = dumE - nu_E_FF;
        membraneI(:,t) = dumI - nu_I_FF;

        
    end

end
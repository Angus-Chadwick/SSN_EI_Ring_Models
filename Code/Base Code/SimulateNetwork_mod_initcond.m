function [rE, rI] = SimulateNetwork_mod_initcond(network, inputs, Nt, NoiseModel,initcond,dt)

NE = network.cells.NE;
NI = network.cells.NI;

rE = zeros([NE, Nt]);
rI = zeros([NI, Nt]);
rE(:,1) = initcond(1:NE);
rI(:,1) = initcond(1:NI);

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

        rE(:,t+1) = rE(:,t) + drE * dt;
        rI(:,t+1) = rI(:,t) + drI * dt;

    end

end
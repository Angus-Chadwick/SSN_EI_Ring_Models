function [rE, rI] = SimulateNetwork(IE_FF, IE_TD, II_FF, II_TD, JEE, JEI, JIE, JII, noise, gamma, tauE, tauI, Nt, NoiseModel)

NE = size(JEE,1);
NI = size(JIE,1);

rE = zeros([NE, Nt]);
rI = zeros([NI, Nt]);

if strcmp(NoiseModel, 'Mult');

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
    noiseI = noiseE/2.5;

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
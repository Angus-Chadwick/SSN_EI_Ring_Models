%% Ring model with attention

clear all

N_TD = 1;

NE = 1000;
NI = NE / 5;


Nloop = 100;


invariant = 'area';

% Time constants

tauE = 10;
tauI = tauE / 2;

Nt = 10000;
dt = tauE / 100;

gamma = 2;

noise = 1;
p=1;



%%

for m=1:3
    
stimvals = 2*pi * [160, 200] / 360;
Nstim = length(stimvals); 
       
RE_TD = zeros([NE, N_TD, Nstim]) ;
RI_TD = zeros([NI, N_TD, Nstim]) ;
RE0 = RE_TD;
RI0 = RI_TD;

RE_TD_std = RE_TD;
RI_TD_std = RI_TD;
RE0_std = RE0;
RI0_std = RI0;


    for q=1:Nstim
        
theta_pE = linspace(0, 2*pi, NE+1);
theta_pE = theta_pE(1:(end-1)); % for circular invariance
theta_pI = linspace(0, 2*pi, NI+1);
theta_pI = theta_pI(1:(end-1));

[~, id1]  = min(abs(stimvals(1) - theta_pI));  % centre stimuli on a cell (centred both on an E and an I due to integer NE/NI)
[~, id2]  = min(abs(stimvals(2) - theta_pI));
stimvals = [theta_pI(id1), theta_pI(id2)];

theta_a = (stimvals(1) + stimvals(2))/2 + pi;

theta_s = stimvals(q);

%% Connectivity

JEE = zeros(NE);
JIE = zeros([NI, NE]);
JEI = zeros([NE, NI]);
JII = zeros(NI);

if m==1 % no tuning

kEE = 1.0;
kIE = kEE *  0; 
kEI = kEE *  0;
kII = +kEE * 0.0 ;

elseif m==2 % iso-orientation inhibition

kEE = 1.0;
kIE = kEE *  0.3; 
kEI = kEE *  0.3;
kII = +kEE * 0.0 ;

elseif m == 3 % cross-orientation inhibition

kEE = 1.0;
kIE = kEE *  0.3;
kEI = -kEE *  0.3;
kII = +kEE * 0.0 ;

end
    
JEE_max = 15/NE;

for i=1:NE

    JEE(i,:) = JEE_max * exp(kEE * cos(theta_pE(i) - theta_pE));
    
end


JEI_mean = mean(JEE(:)) * 2 * 1.1 ;  
JIE_mean = mean(JEE(:)) * 2 * 1.1 ;
JII_mean = mean(JEE(:)) * 1 * 1.1;


JEI_max = JEI_mean / besseli(0, abs(kEI));  % fixed area
JIE_max = JIE_mean / besseli(0, abs(kIE));
JII_max = JII_mean / besseli(0, abs(kII));



for i=1:NE

    JIE(:,i) = JIE_max * exp(kIE * cos(theta_pE(i) - theta_pI));
    JEI(i,:) = JEI_max * exp(kEI * cos(theta_pE(i) - theta_pI));
    
end

for i=1:NI
    
    JII(i,:) = JII_max * exp(kII * cos(theta_pI(i) - theta_pI));

end


%% Stimulus drive

IE_FF_area = 0.005 * 100;  %% seems to implement a gain function
kE_FF = 0.1*kEE;
IE_FF = (IE_FF_area / (2*pi* besseli(0,kE_FF))) * exp(kE_FF * cos(theta_pE - theta_s))';  

II_FF = -0 * ones([NI,1]); %% seems to implement subtractive normalisation

IE_TD_area = IE_FF_area / 10 * 0.5 * p;
kE_TD = 0.05 * 20 * kEE * 0;

II_TD_area = IE_FF_area / 10 * 0.5 * p;
kI_TD = kE_TD * 0.5 * 0;

TD = 'Inh';

if strmatch(TD, 'Exc')

    IE_TD = (IE_TD_area / (2*pi* besseli(0,kE_TD))) * exp(kE_TD * cos(theta_pE - theta_a))';
    II_TD = 0.00000 * ones(size(II_FF));

elseif strmatch(TD, 'Inh')
    
    IE_TD =  0.00000 * ones(size(IE_FF));
    II_TD = II_TD_area / (2*pi* besseli(0,kE_TD)) * exp(kI_TD * cos(theta_pI - theta_a))';

end
 
%% simulate

NoiseModel = 'Add';


parfor n=1:Nloop

    [rE, rI]       = SimulateNetwork(IE_FF, 0*IE_TD, II_FF, 0*II_TD, JEE, JEI, JIE, JII, noise, gamma, tauE, tauI, Nt, NoiseModel);

    
    RE{n}(:,q) = mean(rE(:,300:end),2);
    RI{n}(:,q) = mean(rI(:,300:end),2);
    
    RE_std{n}(:,q) = std(rE(:,300:end),[],2);
    RI_std{n}(:,q) = std(rI(:,300:end),[],2);


    RE_cov{n}(:,:,q) = cov([rE(:,300:end); rI(:,300:end)]')

end

    



RE0 = mean(cat(3,RE{:}),3);
RI0 = mean(cat(3,RI{:}),3);

RE0_std = mean(cat(3,RE_std{:}),3);
RI0_std = mean(cat(3,RI_std{:}),3);

RE0_cov = mean(cat(4,RE_cov{:}),4);


end





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
       
    
end

Xr = rgb2hsv([1,0,0]);
Xb = rgb2hsv([0,0,1]);
Xr(2) = 0.6;
Xb(2) = 0.6;
Xr = hsv2rgb(Xr);
Xb = hsv2rgb(Xb);

Xr2 = rgb2hsv([1,0,0]);
Xb2 = rgb2hsv([0,0,1]);
Xr2(2) = 0.2;
Xb2(2) = 0.2;
Xr2 = hsv2rgb(Xr2);
Xb2 = hsv2rgb(Xb2);



figure

subplot(3,2,1)
hold on
set(gca, 'fontsize', 24)
ylabel('Mean Response')
title('Pyramidal Cells')

plot(theta_pE * 180/pi, RE0_vert{1}, 'color', Xr, 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_vert{1}(1:50:end), 'color', Xr, 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_vert{2}, 'color',Xr2, 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_vert{2}(1:50:end), 'color', Xr2, 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_vert{3}, 'color', 'r', 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_vert{3}(1:50:end), 'color', 'r', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_ang{1}, 'color', Xb, 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_ang{1}(1:50:end), 'color', Xb, 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_ang{2}, 'color', Xb2, 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_ang{2}(1:50:end), 'color', Xb2, 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 2)

plot(theta_pE * 180/pi, RE0_ang{3}, 'color', 'b', 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_ang{3}(1:50:end), 'color', 'b', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)


axis([0,360, 0, 1.1 * max(RE0_vert{3})])
set(gca, 'xtick', 0:90:360)
set(gca, 'xticklabel', {})
box on

subplot(3,2,2)
hold on
set(gca, 'fontsize', 24)
ylabel('Mean Response')
title('Interneurons')


plot(theta_pI * 180/pi, RI0_vert{1}, 'color', Xr, 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_vert{1}(1:10:end), 'color', Xr, 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pI * 180/pi, RI0_vert{2}, 'color',Xr2, 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_vert{2}(1:10:end), 'color', Xr2, 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 3)

plot(theta_pI * 180/pi, RI0_vert{3}, 'color', 'r', 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_vert{3}(1:10:end), 'color', 'r', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)


plot(theta_pI * 180/pi, RI0_ang{1}, 'color', Xb, 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_ang{1}(1:10:end), 'color', Xb, 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pI * 180/pi, RI0_ang{2}, 'color', Xb2, 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_ang{2}(1:10:end), 'color', Xb2, 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 2)

plot(theta_pI * 180/pi, RI0_ang{3}, 'color', 'b', 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_ang{3}(1:10:end), 'color', 'b', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)


axis([0,360, 0, 1.1 * max(RI0_vert{3})])
set(gca, 'xtick', 0:90:360)
set(gca, 'xticklabel', {})
box on

subplot(3,2,3)
hold on
set(gca, 'fontsize', 24)
ylabel('Response std')
title('Pyramidal Cells')


plot(theta_pE * 180/pi, RE0_std_vert{1}, 'color', Xr, 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_std_vert{1}(1:50:end), 'color', Xr, 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_std_vert{2}, 'color',Xr2, 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_std_vert{2}(1:50:end), 'color', Xr2, 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_std_vert{3}, 'color', 'r', 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_std_vert{3}(1:50:end), 'color', 'r', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_std_ang{1}, 'color', Xb, 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_std_ang{1}(1:50:end), 'color', Xb, 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pE * 180/pi, RE0_std_ang{2}, 'color', Xb2, 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_std_ang{2}(1:50:end), 'color', Xb2, 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 2)

plot(theta_pE * 180/pi, RE0_std_ang{3}, 'color', 'b', 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, RE0_std_ang{3}(1:50:end), 'color', 'b', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)



axis([0,360, 0, 1.1 * max(RE0_std_vert{3})])
set(gca, 'xtick', 0:90:360)
set(gca, 'xticklabel', {})
box on

subplot(3,2,4)
hold on
set(gca, 'fontsize', 24)
ylabel('Response std')
title('Interneurons')


plot(theta_pI * 180/pi, RI0_std_vert{1}, 'color', Xr, 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_std_vert{1}(1:10:end), 'color', Xr, 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pI * 180/pi, RI0_std_vert{2}, 'color',Xr2, 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_std_vert{2}(1:10:end), 'color', Xr2, 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 3)

plot(theta_pI * 180/pi, RI0_std_vert{3}, 'color', 'r', 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_std_vert{3}(1:10:end), 'color', 'r', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)


plot(theta_pI * 180/pi, RI0_std_ang{1}, 'color', Xb, 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_std_ang{1}(1:10:end), 'color', Xb, 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pI * 180/pi, RI0_std_ang{2}, 'color', Xb2, 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_std_ang{2}(1:10:end), 'color', Xb2, 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 2)

plot(theta_pI * 180/pi, RI0_std_ang{3}, 'color', 'b', 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, RI0_std_ang{3}(1:10:end), 'color', 'b', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)


axis([0,360, 0, 1.1 * max(RI0_std_vert{3})])
set(gca, 'xtick', 0:90:360)
set(gca, 'xticklabel', {})
box on

subplot(3,2,5)
hold on
set(gca, 'fontsize', 24)
ylabel('Response Selectivity')
title('Pyramidal Cells')

plot(theta_pE * 180/pi, abs(SI0E{1}), 'color', [0.4, 0.4, 0.4], 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, abs(SI0E{1}(1:50:end)), 'color', [0.4,0.4,0.4], 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pE * 180/pi, abs(SI0E{2}), 'color', [0.8, 0.8, 0.8], 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, abs(SI0E{2}(1:50:end)), 'color', [0.8,0.8,0.8], 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 3)

plot(theta_pE * 180/pi, abs(SI0E{3}), 'k', 'linewidth', 3)
plot(theta_pE(1:50:end) * 180/pi, abs(SI0E{3}(1:50:end)), 'color', 'k', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)

axis([0,360, -0.5, 2.5])
set(gca, 'xtick', 0:90:360)
xlabel('Neuron Angle on Ring')
box on

subplot(3,2,6)
hold on
set(gca, 'fontsize', 24)
ylabel('Response Selectivity')
title('Interneurons')


plot(theta_pI * 180/pi, abs(SI0I{1}), 'color', [0.4, 0.4, 0.4], 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, abs(SI0I{1}(1:10:end)), 'color', [0.4,0.4,0.4], 'linestyle', 'none', 'marker', '*', 'markersize', 15, 'linewidth', 3)

plot(theta_pI * 180/pi, abs(SI0I{2}), 'color', [0.8, 0.8, 0.8], 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, abs(SI0I{2}(1:10:end)), 'color', [0.8,0.8,0.8], 'linestyle', 'none', 'marker', '+', 'markersize', 15, 'linewidth', 3)

plot(theta_pI * 180/pi, abs(SI0I{3}), 'k', 'linewidth', 3)
plot(theta_pI(1:10:end) * 180/pi, abs(SI0I{3}(1:10:end)), 'color', 'k', 'linestyle', 'none', 'marker', 'x', 'markersize', 15, 'linewidth', 3)


axis([0,360, -0.5, 2.5])
set(gca, 'xtick', 0:90:360)
box on
xlabel('Neuron Angle on Ring') 

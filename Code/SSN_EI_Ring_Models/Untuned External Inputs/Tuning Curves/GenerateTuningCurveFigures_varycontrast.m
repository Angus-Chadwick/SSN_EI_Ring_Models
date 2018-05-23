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
    
stimvals = pi;
Nstim = 2; % number of stimulus contrasts 
       
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

theta_s = stimvals;
theta_a = 0;

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

IE_FF_area = 0.005 * 100 * q;  %% seems to implement a gain function
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




iv1 = 1; iv2 = 2;

        RE0_std_low{m} = squeeze(RE0_std(:,iv1));
        RE0_std_high{m} = squeeze(RE0_std(:,iv2));

        RE0_low{m} = squeeze(RE0(:,iv1));
        RE0_high{m} = squeeze(RE0(:,iv2));
     
        RI0_std_low{m} = squeeze(RI0_std(:,iv1));
        RI0_std_high{m} = squeeze(RI0_std(:,iv2));

        RI0_low{m} = squeeze(RI0(:,iv1));
        RI0_high{m} = squeeze(RI0(:,iv2));

       

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
subplot(2,3,1)
hold on
set(gca, 'fontsize', 24)
ylabel('Mean Response')
title('Pyramidal Cells')
plot(theta_pE * 180/pi, RE0_low{1}, 'color', Xr, 'linewidth', 4, 'linestyle', '--')
plot(theta_pE * 180/pi, RE0_low{2}, 'color',Xr2, 'linewidth', 4)
plot(theta_pE * 180/pi, RE0_low{3}, 'color', 'r', 'linewidth', 4)
plot(theta_pE * 180/pi, RE0_high{1}, 'color', Xb, 'linewidth', 4, 'linestyle', '--')
plot(theta_pE * 180/pi, RE0_high{2}, 'color', Xb2, 'linewidth', 4)
plot(theta_pE * 180/pi, RE0_high{3}, 'color', 'b', 'linewidth', 4)
set(gca, 'xtick', 0:90:360)
set(gca, 'xticklabel', {})
box on
subplot(2,3,4)
hold on
set(gca, 'fontsize', 24)
ylabel('Mean Response')
title('Interneurons')
plot(theta_pI * 180/pi, RI0_low{1}, 'color', Xr, 'linewidth', 4, 'linestyle', '--')
plot(theta_pI * 180/pi, RI0_low{2}, 'color',Xr2, 'linewidth', 4)
plot(theta_pI * 180/pi, RI0_low{3}, 'color', 'r', 'linewidth', 4)
plot(theta_pI * 180/pi, RI0_high{1}, 'color', Xb, 'linewidth', 4, 'linestyle', '--')
plot(theta_pI * 180/pi, RI0_high{2}, 'color', Xb2, 'linewidth', 4)
plot(theta_pI * 180/pi, RI0_high{3}, 'color', 'b', 'linewidth', 4)
legend('Non-Specific Inhibition', 'Iso-Orientation Inhibition', 'Cross-Orientation Inhibition', '', '', '')
set(gca, 'xtick', 0:90:360)
%set(gca, 'xticklabel', {})
box on
subplot(2,3,2)
hold on
set(gca, 'fontsize', 24)
ylabel('Standard Deviation of Response')
title('Pyramidal Cells')
plot(theta_pE * 180/pi, RE0_std_low{1}, 'color', Xr, 'linewidth', 4, 'linestyle', '--')
plot(theta_pE * 180/pi, RE0_std_low{2}, 'color',Xr2, 'linewidth', 4)
plot(theta_pE * 180/pi, RE0_std_low{3}, 'color', 'r', 'linewidth', 4)
plot(theta_pE * 180/pi, RE0_std_high{1}, 'color', Xb, 'linewidth', 4, 'linestyle', '--')
plot(theta_pE * 180/pi, RE0_std_high{2}, 'color', Xb2, 'linewidth', 4)
plot(theta_pE * 180/pi, RE0_std_high{3}, 'color', 'b', 'linewidth', 4)
set(gca, 'xtick', 0:90:360)
set(gca, 'xticklabel', {})
box on
subplot(2,3,5)
hold on
set(gca, 'fontsize', 24)
ylabel('Standard Deviation of Response')
title('Interneurons')
plot(theta_pI * 180/pi, RI0_std_low{1}, 'color', Xr, 'linewidth', 4, 'linestyle', '--')
plot(theta_pI * 180/pi,  RI0_std_low{2}, 'color',Xr2, 'linewidth', 4)
plot(theta_pI * 180/pi, RI0_std_low{3}, 'color', 'r', 'linewidth', 4)
plot(theta_pI * 180/pi, RI0_std_high{1}, 'color', Xb, 'linewidth', 4, 'linestyle', '--')
plot(theta_pI * 180/pi, RI0_std_high{2}, 'color', Xb2, 'linewidth', 4)
plot(theta_pI * 180/pi, RI0_std_high{3}, 'color', 'b', 'linewidth', 4)
set(gca, 'xtick', 0:90:360)
%set(gca, 'xticklabel', {})
box on
subplot(2,3,3)
hold on
set(gca, 'fontsize', 24)
ylabel('Ratio of Responses')
plot(theta_pE * 180/pi, RE0_high{1} ./ RE0_low{1}, 'color', Xr, 'linewidth', 4, 'linestyle', '--')
plot(theta_pE * 180/pi, RE0_high{2} ./ RE0_low{2}, 'color', Xr2, 'linewidth', 4, 'linestyle', '--')
plot(theta_pE * 180/pi, RE0_high{3} ./ RE0_low{3}, 'color', 'r', 'linewidth', 4, 'linestyle', '--')
subplot(2,3,6)
hold on
set(gca, 'fontsize', 24)
ylabel('Ratio of Responses')
plot(theta_pI * 180/pi, RI0_high{1} ./ RI0_low{1}, 'color', Xr, 'linewidth', 4, 'linestyle', '--')
plot(theta_pI * 180/pi, RI0_high{2} ./ RI0_low{2}, 'color', Xr2, 'linewidth', 4, 'linestyle', '--')
plot(theta_pI * 180/pi, RI0_high{3} ./ RI0_low{3}, 'color', 'r', 'linewidth', 4, 'linestyle', '--')

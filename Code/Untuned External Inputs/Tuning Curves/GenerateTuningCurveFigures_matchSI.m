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

if m==1

JEI_mean = 0.04 ;  
JIE_mean = 0.04 ;
JII_mean = mean(JEE(:)) * 1 * 1.1;

elseif m==2
    
JEI_mean = 0.035 ;  
JIE_mean = 0.035 ;
JII_mean = mean(JEE(:)) * 1 * 1.1;

elseif m==3 
    
JEI_mean = 0.045 ;  
JIE_mean = 0.045 ;
JII_mean = mean(JEE(:)) * 1 * 1.1;

end

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
        
        ExternalE_Input{m} = IE_FF;
        RecurrentEE_Input{m} = JEE * RE0_ang{m};
        RecurrentEI_Input{m} = JEI * RI0_ang{m};
    
        ExternalI_Input{m} = II_FF;
        RecurrentIE_Input{m} = JIE * RE0_ang{m};
        RecurrentII_Input{m} = JII * RI0_ang{m};
    
end

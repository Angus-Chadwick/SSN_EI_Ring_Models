%% Ring model with attention

clear all

N_TD = 1;

NE = 1000;
NI = NE / 5;


Nloop = 100;

RE_covtot = cell(Nloop);

invariant = 'area';

% Time constants

tauE = 10;
tauI = tauE / 2;

Nt = 10000;
dt = tauE / 100;

gamma = 2;

noise = 1;


stimvals = 2*pi * [160,200] / 360;
Nstim = length(stimvals); 
       


RE = cell([Nloop, 1]);
RI = RE;

RE_std = RE;
RI_std = RI;
RE_cov = RE;
Rtot_cov = RE;

p = 360;

for q=1:Nstim
        
theta_pE = linspace(0, 2*pi, NE+1);
theta_pE = theta_pE(1:(end-1)); % for circular invariance
theta_pI = linspace(0, 2*pi, NI+1);
theta_pI = theta_pI(1:(end-1));

[~, id1]  = min(abs(stimvals(1) - theta_pI));  % centre stimuli on a cell (centred both on an E and an I due to integer NE/NI)
[~, id2]  = min(abs(stimvals(2) - theta_pI));
stimvals = [theta_pI(id1), theta_pI(id2)];

theta_a = (stimvals(1) + stimvals(2))/2;

theta_s = stimvals(q);

%% Connectivity

JEE = zeros(NE);
JIE = zeros([NI, NE]);
JEI = zeros([NE, NI]);
JII = zeros(NI);

kEE = 1.0;
kIE = 0.5; % 20, 101
kEI = -0.5;
kII = +kEE * 0.0 ;



JEE_max = 15/NE;

for i=1:NE

    JEE(i,:) = JEE_max * exp(kEE * cos(theta_pE(i) - theta_pE));
    
end


JEI_mean = mean(JEE(:)) * 2 * 1.1 ;  
JIE_mean = mean(JEE(:)) * 2 * 1.1 ;
JII_mean = mean(JEE(:)) * 1 * 1.1;


if strmatch(invariant, 'area')

JEI_max = JEI_mean / besseli(0, abs(kEI));  % fixed area
JIE_max = JIE_mean / besseli(0, abs(kIE));
JII_max = JII_mean / besseli(0, abs(kII));

elseif strmatch(invariant, 'maxval')

JEI_max = JEI_mean /  besseli(0, 0)  * (exp(0.5) / exp(abs(kEI)));   % fixed maximum
JIE_max = JIE_mean / besseli(0, 0)  * (exp(0.5) / exp(abs(kIE)));  
JII_max = JII_mean / besseli(0, 0)  * (exp(0.5) / exp(abs(kII)));

end


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

IE_TD_area = IE_FF_area / 10 * 0.5 * (p-1);
kE_TD = 0.05 * 20 * kEE * 0;

II_TD_area = IE_FF_area / 100 * (p-400);
kI_TD = 0.5;

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
    [rE_TD, rI_TD]       = SimulateNetwork(IE_FF, IE_TD, II_FF, II_TD, JEE, JEI, JIE, JII, noise, gamma, tauE, tauI, Nt, NoiseModel);

    
    RE{n}(:,q) = mean(rE(:,300:end),2);
    RI{n}(:,q) = mean(rI(:,300:end),2);

    RE_TD{n}(:,q) = mean(rE_TD(:,300:end),2);
    RI_TD{n}(:,q) = mean(rI_TD(:,300:end),2);
    
    RE_std{n}(:,q) = std(rE(:,300:end),[],2);
    RI_std{n}(:,q) = std(rI(:,300:end),[],2);

    RE_std_TD{n}(:,q) = std(rE_TD(:,300:end),[],2);
    RI_std_TD{n}(:,q) = std(rI_TD(:,300:end),[],2);
    
    RE_cov{n}(:,:,q) = cov(rE(:,300:end)');
    Rtot_cov{n}(:,:,q) = cov([rE(:,300:end); rI(:,300:end)]')
    
    RE_cov_TD{n}(:,:,q) = cov(rE_TD(:,300:end)');
    Rtot_cov_TD{n}(:,:,q) = cov([rE_TD(:,300:end); rI(:,300:end)]')
    
end

end

RE0 = mean(cat(3,RE{:}),3);
RI0 = mean(cat(3,RI{:}),3);

RE0_std = mean(cat(3,RE_std{:}),3);
RI0_std = mean(cat(3,RI_std{:}),3);

RE0_cov = mean(cat(4,RE_cov{:}),4);
Rtot0_cov = mean(cat(4,Rtot_cov{:}),4);

RE0_TD = mean(cat(3,RE_TD{:}),3);
RI0_TD = mean(cat(3,RI_TD{:}),3);

RE0_std_TD = mean(cat(3,RE_std_TD{:}),3);
RI0_std_TD = mean(cat(3,RI_std_TD{:}),3);

RE0_cov_TD = mean(cat(4,RE_cov_TD{:}),4);
Rtot0_cov_TD = mean(cat(4,Rtot_cov_TD{:}),4);




mmE = zeros(Nstim);
mmI = zeros(Nstim);
m0E = zeros(Nstim);
m0I = zeros(Nstim);

iv1 = 1;iv2 = 2;

        RE0_std_vert = squeeze(RE0_std(:,iv1));
        RE0_std_ang = squeeze(RE0_std(:,iv2));

        RE0_vert = squeeze(RE0(:,iv1));
        RE0_ang = squeeze(RE0(:,iv2));
     
        RI0_std_vert = squeeze(RI0_std(:,iv1));
        RI0_std_ang = squeeze(RI0_std(:,iv2));

        RI0_vert = squeeze(RI0(:,iv1));
        RI0_ang  = squeeze(RI0(:,iv2));

        % E selectivity
        
        dmu_0E = RE0_ang - RE0_vert;
        
        poolvar_0E = (0.5* (RE0_std_ang.^2 + RE0_std_vert.^2));

        SI0E = dmu_0E.^2 ./ poolvar_0E;

        m0E(iv1,iv2) = nansum(SI0E);

        % I selectivity
        
        dmu_0I = RI0_ang - RI0_vert;
        
        poolvar_0I = (0.5* (RI0_std_ang.^2 + RI0_std_vert.^2));

        SI0I = dmu_0I.^2 ./ poolvar_0I;
        
        % TD
        
        RE0_std_vert_TD = squeeze(RE0_std_TD(:,iv1));
        RE0_std_ang_TD = squeeze(RE0_std_TD(:,iv2));

        RE0_vert_TD = squeeze(RE0_TD(:,iv1));
        RE0_ang_TD = squeeze(RE0_TD(:,iv2));
     
        RI0_std_vert_TD = squeeze(RI0_std_TD(:,iv1));
        RI0_std_ang_TD = squeeze(RI0_std_TD(:,iv2));

        RI0_vert_TD = squeeze(RI0_TD(:,iv1));
        RI0_ang_TD  = squeeze(RI0_TD(:,iv2));

        % E selectivity
        
        dmu_0E_TD = RE0_ang_TD - RE0_vert_TD;
        
        poolvar_0E_TD = (0.5* (RE0_std_ang_TD.^2 + RE0_std_vert_TD.^2));

        SI0E_TD = dmu_0E_TD.^2 ./ poolvar_0E_TD;



RE_covtot1{p} = RE0_cov(:,:,1);
RE_covtot2{p} = RE0_cov(:,:,2);
SItot_E(p) = squeeze(RE0(:,1) - RE0(:,2))' * inv(0.5 * (RE_covtot1{p} + RE_covtot2{p})) * squeeze(RE0(:,1) - RE0(:,2)); 
SItot_E_ind(p) =  squeeze(RE0(:,1) - RE0(:,2))' * inv(0.5 * diag(diag(RE_covtot1{p} + RE_covtot2{p}))) * squeeze(RE0(:,1) - RE0(:,2)); 
SI_I(p) = m0I(1,2);
SI_E(p) = m0E(1,2);

RE_covtot1_TD{p} = RE0_cov_TD(:,:,1);
RE_covtot2_TD{p} = RE0_cov_TD(:,:,2);
SItot_E_TD(p) = squeeze(RE0_TD(:,1) - RE0_TD(:,2))' * inv(0.5 * (RE_covtot1_TD{p} + RE_covtot2_TD{p})) * squeeze(RE0_TD(:,1) - RE0_TD(:,2)); 
SItot_E_ind_TD(p) =  squeeze(RE0_TD(:,1) - RE0_TD(:,2))' * inv(0.5 * diag(diag(RE_covtot1_TD{p} + RE_covtot2_TD{p}))) * squeeze(RE0_TD(:,1) - RE0_TD(:,2)); 


subplot(2,2,3)
hold on
plot(RE0(:,1), 'color', 'r', 'linewidth', 4)
plot(RE0(:,2), 'color', 'b', 'linewidth', 4)
plot(RE0_TD(:,1), 'color', 'r', 'linewidth', 4)
plot(RE0_TD(:,2), 'color', 'b', 'linewidth', 4)
set(gca,'fontsize', 24)
set(gca, 'xtick', [])
xlabel('Neuron Location on Ring')
ylabel('Response Amplitude')
box on
subplot(2,2,4)
hold on
plot(sqrt(SI0E), 'color', 'k', 'linewidth', 4)
plot(sqrt(SI0E_TD), 'color', 'k', 'linewidth', 4)
set(gca,'fontsize', 24)
set(gca, 'xtick', [])
xlabel('Neuron Location on Ring')
ylabel('Selectivity')
box on

subplot(2,1,1)

% plot SIind and SIpop vs top-down input here

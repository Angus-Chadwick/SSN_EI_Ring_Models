%% This script is used to find fixed points of the Supralinear Stabilised Network with arbitrary tuning of all 4 connection types and inputs

clear all

% simulation parameters

Nvals = 11; % 
Nt = 1000;  % number of timesteps

% initialise variables

RE_covtot = cell(Nvals);
RE = cell(1);
RI = RE;
RE_std = RE;
RI_std = RI;
RE_cov = RE;
Rtot_cov = RE;

% fixed input parameters
NE = 1000;

stimvals = 2*pi * [180,220] / 360;
Nstim = length(stimvals);

noise = 0;
kE_FF = 0.5; % standard value = 0.5
kI_FF = 0.0;
IE_FF_area = 0.005 * 100;
JEE_mean = 15/NE * besseli(0,1);
JEI_mean = 0.04;  % standard values
JIE_mean = 0.04;
JII_mean = JEE_mean * 1.1;

stepE = 25 /5;
stepI = 50 /5;

ParameterSweep = 'NonUniformRing';
soft_thresh = 0;
for nEI = 1:(2*Nvals)
    nEI
    for nIE = 1:(2*Nvals)
   
q=1;

%% create network        

if strcmp(ParameterSweep, 'NoisyWeights')

kEE = 3;
kIE = 0.1 * (nEI-1); 
kEI = 0.4; % I to E concentration
kII = 0.0; 

II_FF_area = IE_FF_area * (0);

noiseEE = 0;
noiseIE = (nIE-1) / 10;
noiseEI = 0;
noiseII = 0;

elseif strcmp(ParameterSweep, 'NonUniformRing')

kEE = 2;  % note: now trying this reduced in all sims... (used to be 2)
kIE = 0.1; 
kEI = 0.4; % I to E concentration
kII = 0.0; 

II_FF_area = IE_FF_area * (0);

elseif strcmp(ParameterSweep, 'EvsIWeights')

kEE = (nEI-1) / stepE * 2.5;   % for sweepss
kIE = (nIE-1) / stepI * 2;  % for sweeps
% kEE = (nEI-1) / stepE; 
% kIE = (nIE-1) / stepI;
kEI = 0.4; % I to E concentration
kII = 0.0; 

II_FF_area = IE_FF_area * (0);

elseif strcmp(ParameterSweep, 'EvsIAmps')

kEE = 2;
kIE = 0.1;
kEI = 0.4;
kII = 0.0;

II_FF_area = IE_FF_area * (0);

JEE_mean = 0.04  * nEI / 22;
JIE_mean = 0.08  * nIE / 22;

elseif strcmp(ParameterSweep, 'EEvsIEWeights')

kEE = kE_FF * (nIE-1) / stepE;
kIE = kEE; 
kEI = kE_FF * (nEI-1) / stepE; % I to E concentration
kII = 3 * kEI; 

II_FF_area = 0;

JEE_mean = 0.02;
JEI_mean = 0.015;
JIE_mean = 0.015;
JII_mean = 0.03;

elseif strcmp(ParameterSweep, 'WeightsConcvsAmp')

    kE_FF = 2;
    kI_FF = 2;
    
kEE = kE_FF * 2;
kIE = kEE; 
kEI = kEE; % I to E concentration
kII = kEE; 

II_FF_area = - IE_FF_area;
JEE_mean = 15/NE * besseli(0,1);
JEI_mean = 0.024;  % standard values
JIE_mean = 0.024 * (1 - nEI / (2*Nvals+1));
JII_mean = JEE_mean * 1.1;

% 
% JEE_max = 0.02 * nEI / stepI;
% JEI_max = 0.04 * nEI / stepI;  % standard values
% JIE_max = 0.04 * nEI / stepI;
% JII_max = 0.02 * nEI / stepI;
% 
% JEE_mean = JEE_max * besseli(0, abs(kEE)) ;
% JIE_mean = JIE_max * besseli(0, abs(kIE));
% JEI_mean = JEI_max * besseli(0, abs(kEI));
% JII_mean = JII_max * besseli(0, abs(kII));

elseif strcmp(ParameterSweep, 'WeightsAmp')
% 
%    % broader I
%     
% kEE = kE_FF * 3;
% kIE = kE_FF * 2; 
% kEI = kE_FF * 3; % I to E concentration
% kII = kE_FF * 3; 
% 
% II_FF_area = 0;
% 
% JEE_mean = 0.0165;
% JIE_mean = 0.02 ;
% JEI_mean = 0.02 ;
% JII_mean = 0.0165 ;
% 
%    % narrower I
%     
% kEE = kE_FF * 3;
% kIE = kE_FF * 3; 
% kEI = kE_FF * 3; % I to E concentration
% kII = kE_FF * 3; 
% 
% II_FF_area = 0;
% 
% JEE_mean = 0.0165;
% JIE_mean = 0.02 ;
% JEI_mean = 0.02 ;
% JII_mean = 0.0165 ;
% 
%    % broader I
%     
% kEE = 3;
% kIE = kE_FF; 
% kEI = kE_FF; % I to E concentration
% kII = kE_FF; 
% 
% II_FF_area = 0;
% 
% JEE_mean = 0.019;
% JIE_mean = 0.04 ;
% JEI_mean = 0.04 ;
% JII_mean = 0.019 ;
% 
% 
%    % narrower I
%     
% kEE = 3;
% kIE = kE_FF * 2; 
% kEI = kE_FF * 2; % I to E concentration
% kII = kE_FF * 2; ca
% 
% II_FF_area = 0;
% 
% JEE_mean = 0.019;
% JIE_mean = 0.0325 ;
% JEI_mean = 0.0325 ;
% JII_mean = 0.019 ;
% 
%    % even narrower I  
% 
% kEE = 3;
% kIE = kE_FF * 4; 
% kEI = kE_FF * 4; % I to E concentration
% kII = kE_FF * 4; 
% 
% II_FF_area = 0;
% 
% JEE_mean = 0.019;
% JIE_mean = 0.0245 ;
% JEI_mean = 0.0245;
% JII_mean = 0.019 ;
% 
% 
% 
%    % scale whole system
% 
% kEE = 3/2;
% kIE = kE_FF * 4 / 2; 
% kEI = kE_FF * 4 / 2; % I to E concentration
% kII = kE_FF * 4 / 2; 
% 
% II_FF_area = 0;
% 
% JEE_mean = 0.019;
% JIE_mean = 0.02375 ;
% JEI_mean = 0.02375;
% JII_mean = 0.019 ;
% 
%    % scale whole system further
% 
% kEE = 1;
% kIE = kE_FF * 4 / 3; 
% kEI = kE_FF * 4 / 3; % I to E concentration
% kII = kE_FF * 4 / 3; 
% 
% II_FF_area = 0;
% 
% JEE_mean = 0.019;
% JIE_mean = 0.023 ;
% JEI_mean = 0.023;
% JII_mean = 0.019 ;

% broad I (narrower FF than rec weights)

kE_FF = 2;

kEE = 1 ;
kIE = 0.1 ; 
kEI = 0.1 ; % I to E concentration
kII = 0 ; 

II_FF_area = -IE_FF_area;
kI_FF = kE_FF;

JEE_mean = 0.019;
JIE_mean = 0.0285 ;
JEI_mean = 0.0285 ;
JII_mean = 0.02 ;

% narrower I 

kE_FF = 2;

kEE = 1 ;
kIE = 0.3 ; 
kEI = 0.3 ; % I to E concentration
kII = 0 ; 

II_FF_area = -IE_FF_area;
kI_FF = kE_FF;

JEE_mean = 0.019;
JIE_mean = 0.0256 ;
JEI_mean = 0.0256 ;
JII_mean = 0.02 ;

% narrowest

kE_FF = 2;

kEE = 1 ;
kIE = 1 ; 
kEI = 1 ; % I to E concentration
kII = 1 ; 

II_FF_area = -IE_FF_area;
kI_FF = kE_FF;

JEE_mean = 0.019;
JIE_mean = 0.0207 ;
JEI_mean = 0.0207 ;
JII_mean = 0.02 ;

% randomly identified 1

IE_FF_area = 0.5;
II_FF_area = 0.0869;
kE_FF = 2;
kI_FF = 2;

kEE = 3.2542;
kIE = 6.2309;
kEI = 0.1750;
kII = 0.2414;

JEE_mean = 0.0180;
JEI_mean = 0.0275;
JIE_mean = 0.0255;
JII_mean = 0.0190;
 
% randomly identified 2

IE_FF_area = 0.5;
II_FF_area = -0.3207;
kE_FF = 2;
kI_FF = 2;

kEE = 3.8961;
kIE = 1.6084;
kEI = 2.7213;
kII = 3.4284;

JEE_mean = 0.0102;
JEI_mean = 0.1337;
JIE_mean = 0.0055;
JII_mean = 0.0256;




% original paper modelling

% iso

kE_FF = 0.5;
kI_FF = 0;
IE_FF_area = 0.5;
II_FF_area = 0;

kEE = 1;
kIE = 0.5;
kEI = 0.5;
kII = 0;

JEE_mean = 15/NE * besseli(0,1);
JIE_mean = 0.024;
JEI_mean = 0.024;
JII_mean = JEE_mean * 1.1;

% untuned

kE_FF = 0.5;
kI_FF = 0.0;
IE_FF_area = 0.5;
II_FF_area = 0;

kEE = 1;
kIE = 0.1;
kEI = 0.5;
kII = 0;

JEE_mean = 15/NE * besseli(0,1);
JIE_mean = 0.030;
JEI_mean = 0.030;
JII_mean = JEE_mean * 1.1;

elseif strcmp(ParameterSweep, 'random')
E = 0;

kE_FF = 2;
kI_FF = kE_FF;

while E < 2

kEE = kE_FF * max(0, randn+1);
kIE = kE_FF * max(0, randn+1);
kEI = kE_FF * max(0, randn+1); % I to E concentration
kII = kE_FF * max(0, randn+1); 

II_FF_area = randn * IE_FF_area;

JEE_mean = 0.02 * max(0, randn+1);
JIE_mean = 0.04 * max(0, randn+1);
JEI_mean = 0.04 * max(0, randn+1);
JII_mean = 0.02 * max(0, randn+1);

Jlin = [JEE_mean, -JEI_mean; JIE_mean, -JII_mean];
E = sum(real(eig(Jlin)) < 0);

end

randpars.concs.kEE(nEI,nIE) = kEE;
randpars.concs.kIE(nEI,nIE) = kIE;
randpars.concs.kEI(nEI,nIE) = kEI;
randpars.concs.kII(nEI,nIE) = kII;
randpars.amps.JEE_mean(nEI,nIE) = JEE_mean;
randpars.amps.JIE_mean(nEI,nIE) = JIE_mean;
randpars.amps.JEI_mean(nEI,nIE) = JEI_mean;
randpars.amps.JII_mean(nEI,nIE) = JII_mean;
randpars.inputs.II_FF_area(nEI,nIE) = II_FF_area;

elseif strcmp(ParameterSweep, 'Persi')

    kE_FF = 4;
    kI_FF = 1 * kE_FF;
    
    kEE = 2*kE_FF;
    kIE = (kE_FF * kI_FF) / (kE_FF - 0.5 * kI_FF) ;
    kEI = (kE_FF * kI_FF) / (kI_FF - 0.5 * kE_FF) ;
    kII = 2*kI_FF;
    
    IE_FF_area = 10;
    II_FF_area = -IE_FF_area / 10;
    
    JEE_mean = 0.02;
    JIE_mean = 0.024;
    JEI_mean = 0.024;
    JII_mean = 0.02 ;

end

network = create_network_varyall(kEE,kEI,kIE,kII, JEE_mean, JEI_mean, JIE_mean, JII_mean);
NE = network.cells.NE;
NI = network.cells.NI;
if strmatch(ParameterSweep,'NoisyWeights')
network.connectivity.JEE = max(network.connectivity.JEE + mean(network.connectivity.JEE(:)) * randn(NE) * noiseEE, 0);
network.connectivity.JEI = max(network.connectivity.JEI + mean(network.connectivity.JEI(:)) * randn(NE,NI) * noiseEI, 0);
network.connectivity.JIE = max(network.connectivity.JIE + mean(network.connectivity.JIE(:)) * randn(NI,NE) * noiseIE, 0);
network.connectivity.JII = max(network.connectivity.JII + mean(network.connectivity.JII(:)) * randn(NI) * noiseII, 0);
end

%% create inputs

theta_s = stimvals(q);
theta_aI = theta_s;

inputs  = create_inputs_varyall(theta_s, noise, kE_FF, kI_FF, IE_FF_area, II_FF_area, network);

if strmatch(ParameterSweep, 'NonUniformRing')
    
%     kmod = 0.1 * (nIE - 1);
%     Amod = 0.05 * (nEI - 1); %  / besseli(0,kmod);  % trying out bessel normalisation - seems better without...
%     JIE_mod = exp(kmod * cos(inputs.theta_pI - stimvals(1)))' *  exp(kmod * cos(inputs.theta_pE - stimvals(1))) + exp(kmod * cos(inputs.theta_pI - stimvals(2)))' *  exp(kmod * cos(inputs.theta_pE - stimvals(2)));  
%     network.connectivity.JIE = network.connectivity.JIE .* (1+Amod * JIE_mod) ./ mean(1+ Amod * JIE_mod(:));
% 
%     kmod = 0.4 * (nIE - 1);
%     Amod = 0.000075 * (nEI - 1)  / besseli(0,kmod)^2;
%     
%     JIE_mod = exp(kmod * cos(inputs.theta_pI - stimvals(1)))' *  exp(kmod * cos(inputs.theta_pE - stimvals(1))) + exp(kmod * cos(inputs.theta_pI - stimvals(2)))' *  exp(kmod * cos(inputs.theta_pE - stimvals(2)));
%     network.connectivity.JIE = (network.connectivity.JIE + Amod * JIE_mod) .* mean(network.connectivity.JIE(:)) ./ mean(network.connectivity.JIE(:) + Amod * JIE_mod(:));


    kmod = 0.2 * (nIE - 1);
    Amod = 3 * 0.000075 * (nEI - 1)  / besseli(0,kmod)^2; % just one stim because its modes around a fixed point
    
    Nsubnets = 1; % number of stimulus-specific subnetworks
    
    if Nsubnets == 1
    
    JIE_mod = exp(kmod * cos(inputs.theta_pI - stimvals(1)))' *  exp(kmod * cos(inputs.theta_pE - stimvals(1)));

    elseif Nsubnets == 2
    
    JIE_mod = exp(kmod * cos(inputs.theta_pI - stimvals(1)))' *  exp(kmod * cos(inputs.theta_pE - stimvals(1))) + exp(kmod * cos(inputs.theta_pI - stimvals(2)))' *  exp(kmod * cos(inputs.theta_pE - stimvals(2)));
            
    end
        
    network.connectivity.JIE = (network.connectivity.JIE + Amod * JIE_mod) .* mean(network.connectivity.JIE(:)) ./ mean(network.connectivity.JIE(:) + Amod * JIE_mod(:));



%     kmod = 0.1 * (nIE - 1);
%     Amod = 0.0002 * (nEI - 1)  / besseli(0,kmod)^2;
% 
%     JEE_mod = exp(kmod * cos(inputs.theta_pE - stimvals(1)))' *  exp(kmod * cos(inputs.theta_pE - stimvals(1))) + exp(kmod * cos(inputs.theta_pE - stimvals(2)))' *  exp(kmod * cos(inputs.theta_pE - stimvals(2)));
%     network.connectivity.JEE = (network.connectivity.JEE + Amod * JEE_mod) .* mean(network.connectivity.JEE(:)) ./ mean(network.connectivity.JEE(:) + Amod * JEE_mod(:));
%     
end   
%% simulate
 
 NoiseModel = 'Add'; 

[rE, rI]       = SimulateNetwork_mod_varyall(network, inputs, Nt, NoiseModel);

StableSim(nEI,nIE) = ~isnan(sum(rE(:)));

if StableSim(nEI,nIE) R0 = [mean(rE(:,(Nt/2):end),2)', mean(rI(:,(Nt/2):end),2)']; elseif or(strcmp(ParameterSweep, 'random'), nIE == 1) R0 = [inputs.IE_FF', inputs.II_FF']; else R0 = FixedPoint{nEI,nIE-1};end  % set initial guess
    
            
            FixedPointFinder;
            R0 = max(0,rmin');
            Phip = diag(2 * R0.^(1/2));
        
            JacobianType = 'dual';
           
            if strcmp(JacobianType, 'normal')
            
            Jstar = inv(T) * (Phip * W - eye(NE+NI)); % jacobian of output system, inv(T) is in front
        
            elseif strcmp(JacobianType, 'dual') 
            
            Jstar = (W * Phip - eye(NE+NI)) * inv(T); % jacobian of input system, inv(T) is at the back
                
            end
                
            DerivativeNorm(nEI,nIE) = norm(inv(T) * ( R0' - max(0, W * R0' + [inputs.IE_FF', inputs.II_FF']').^2));
            FixedPoint{nEI,nIE} = R0;
        

        
        [Vleft,D] = eig(Jstar');  % try computing left eigenvectors and eigenvalues simply by taking transpose of matrix, no need to invert
        Vleft = Vleft';
        Evals{nEI,nIE} = diag(D);        
        
        inputs.noise = 2;
        Inp = ([inputs.IE_FF .* (- kE_FF * sin(inputs.theta_pE - theta_s))'; inputs.II_FF .* (- kI_FF * sin(inputs.theta_pI - theta_s))']);

        optvec_method = 'SNRdual';
        
        if strcmp(optvec_method, 'SNRdual')
        
        CovInp =  diag([inputs.noise * mean(inputs.IE_FF) * ones(NE,1); inputs.noise/2 * mean(inputs.IE_FF) * ones(NI,1)]); % removed inv(T) from this expression as it is in input basis
            
        optvec =  Inp ./ diag(CovInp); % works only for diagonal input covariance
        optvec(abs(Inp) < 1e-10) = 0;
        optvec = optvec / norm(optvec);        
        
        elseif strcmp(optvec_method, 'SNRnormal')
        
                Inp = Phip * ([inputs.IE_FF .* (- kE_FF * sin(inputs.theta_pE - theta_s))'; zeros(NI,1)]);
                CovInp = inv(T) * Phip *  diag([inputs.noise * mean(inputs.IE_FF) * ones(NE,1); inputs.noise/2 * mean(inputs.IE_FF) * ones(NI,1)]) * inv(T) * Phip;

               optvec =  Inp ./ diag(CovInp); % works only for diagonal input covariance
               optvec(abs(Inp) < 1e-10) = 0;
               optvec = optvec / norm(optvec);    
            
        elseif strcmp(optvec_method, 'Signal')           

            optvec = Inp / norm(Inp);
            
        end
         
        
        p = find(max(real(diag(D))),1);
        SNRmode{nEI,nIE} = (Vleft * Inp).^2 ./ diag(Vleft * CovInp * Vleft') / (Inp' * pinv(CovInp) * Inp);
        taumode{nEI,nIE} = -1./real(diag(D));        
        
        SNRFrac(nEI,nIE) = sqrt(max(real(SNRmode{nEI,nIE})));  % complex translation mode is problematic
        tau(nEI,nIE) = taumode{nEI,nIE}(find(real(SNRmode{nEI,nIE}) == max(real(SNRmode{nEI,nIE})),1));
       
        
        Modes{nEI,nIE} = Vleft;
        
        
        % compute stationary state information
%         
%         J = (Phip * W - eye(NE+NI));
%         Sigma = lyap((inv(T) * J), inv(T) * Phip * CovInp *inv(T) * Phip); 
%         rp = -inv(J) * Phip * Inp;
%         InfOut(nEI,nIE) =  rp' * pinv(Sigma) * rp  / (Inp' * pinv(CovInp) * Inp)

% % SNR for complex modes:
% 
% S =  ((real(D) * real(Vleft) + imag(D) * imag(Vleft)) * Inp) .*  (-2./diag(abs(D)).^2); % signal for oscillating mode (projection onto real part)
% N =  diag(real(Vleft) * CovInp * real(Vleft')) .* (-real(diag(D)) ./ abs(diag(D)).^2 - 1./ real(diag(D))) + diag(imag(Vleft) * CovInp * imag(Vleft')) .* (+real(diag(D)) ./ abs(diag(D)).^2 - 1./ real(diag(D))) + diag(real(Vleft) * CovInp * imag(Vleft')) .* (-2 * imag(diag(D)) ./ abs(diag(D)).^2);
% SNRmode_full{nEI,nIE} = 0.5 * S.^2 ./ N / (Inp' * pinv(CovInp) * Inp);   
% InputSNRs = SNRmode_full{nEI,nIE} ./ taumode{nEI,nIE};
% SNRFrac_full(nEI,nIE) = max(InputSNRs);
% q = find(InputSNRs == max(InputSNRs),1);
% tau_full(nEI,nIE) = taumode{nEI,nIE}(q);

    end
end



%% For image maps

for nIE = 1:22
    for nEI = 1:22
        StableFP(nEI,nIE) = and(StableSim(nEI,nIE), min(taumode{nEI,nIE}) > 0);
        Width(nEI,nIE) = sum(FixedPoint{nEI,nIE}(1:1000) > 0);
    end
end

Sharpening = 1000./Width - 1;

% create a colormap for scatter plots

clear colorMap

for i=1:255
    
            colorMap(i,:) = hsv2rgb([1,1,(i-1)/254]);

end


if strcmp(ParameterSweep, 'EvsIWeights')


y = ([1:22]-1) / stepE * 2.5; 
x = ([1:22]-1) / stepI * 2;

subplot(3,4,1)
Sharpening(StableFP == 0) = nan;
h = imagesc(x, y, Sharpening)
h.AlphaData = ones(size(h.CData)); 
h.AlphaData(isnan(h.CData)) = 0;
set(gca, 'fontsize', 18)
xlabel('E to I Specificity')
ylabel('E to E Specificity')
title('Recurrent Sharpening')
set(gca, 'ydir', 'normal')
colormap(gca, colorMap)
h = colorbar
caxis([0,5])


subplot(3,4,2)
SNRFrac(StableFP == 0) = nan;
h = imagesc(x, y, SNRFrac);
h.AlphaData = ones(size(h.CData)); 
h.AlphaData(isnan(h.CData)) = 0;
set(gca, 'fontsize', 18)
xlabel('E to I Specificity')
ylabel('E to E Specificity')
title('Norm. Mode Input SNR')
h = colorbar;
caxis([0,1])
set(gca, 'ydir', 'normal')

subplot(3,4,3)
tau(StableFP == 0) = nan;
h = imagesc(x, y, tau);
h.AlphaData = ones(size(h.CData)); 
h.AlphaData(isnan(h.CData)) = 0;
set(gca, 'fontsize', 18)
xlabel('E to I Specificity')
ylabel('E to E Specificity')
title('Mode Time Constant')
h = colorbar;
caxis([0,300])
set(gca, 'ydir', 'normal')

subplot(3,4,4)
hold on
for i=1:size(SNRmode,1)
for j=1:size(SNRmode,2)
if StableFP(i,j)
    clr = hsv2rgb([1,1,min(Sharpening(i,j) / 3,1)]);
    scatter(real(sqrt(real(SNRmode{i,j}))), taumode{i,j}, 100, 'marker', '.', 'markeredgecolor', clr, 'linewidth', 5)
end
end
end

set(gca, 'fontsize', 18)
xlabel('Norm. Mode Input SNR')
ylabel('Mode Time Constant')
axis([0,1,0,300])

colormap(gca, colorMap)
h = colorbar
set(h,'fontsize', 18)
title(h,'Sharpening')
set(h, 'ytick', [0,1])
set(h, 'yticklabel', [{'0'}, {'3'}])


elseif strcmp(ParameterSweep, 'EvsIAmps')

y = 0.04  * [1:22] / 22;
x = 0.08  * [1:22] / 22;

subplot(3,4,5)
Sharpening(StableFP == 0) = nan;
h = imagesc(x, y, Sharpening)
h.AlphaData = ones(size(h.CData)); 
h.AlphaData(isnan(h.CData)) = 0;
set(gca, 'fontsize', 18)
xlabel('E to I Weight')
ylabel('E to E Weight')
title('Recurrent Sharpening')
set(gca, 'ydir', 'normal')
colormap(gca, colorMap)
h = colorbar
caxis([0,5])


subplot(3,4,6)
SNRFrac(StableFP == 0) = nan;
h = imagesc(x, y, SNRFrac);
h.AlphaData = ones(size(h.CData)); 
h.AlphaData(isnan(h.CData)) = 0;
set(gca, 'fontsize', 18)
xlabel('E to I Weight')
ylabel('E to E Weight')
title('Norm. Mode Input SNR')
h = colorbar;
caxis([0,1])
set(gca, 'ydir', 'normal')
subplot(3,4,7)
tau(StableFP == 0) = nan;
h = imagesc(x, y, tau);
h.AlphaData = ones(size(h.CData)); 
h.AlphaData(isnan(h.CData)) = 0;
set(gca, 'fontsize', 18)
xlabel('E to I Weight')
ylabel('E to E Weight')
title('Mode Time Constant')
h = colorbar;
caxis([0,300])
set(gca, 'ydir', 'normal')

subplot(3,4,8)
hold on
for i=1:size(SNRmode,1)
for j=1:size(SNRmode,2)
if StableFP(i,j)
    clr = hsv2rgb([1,1,min(Sharpening(i,j) / 3,1)]);
    scatter(real(sqrt(real(SNRmode{i,j}))), taumode{i,j}, 100, 'marker', '.', 'markeredgecolor', clr, 'linewidth', 5)
end
end
end

set(gca, 'fontsize', 18)
xlabel('Norm. Mode Input SNR')
ylabel('Mode Time Constant')
axis([0,1,0,300])

colormap(gca, colorMap)
h = colorbar
set(h,'fontsize', 18)
title(h,'Sharpening')
set(h, 'ytick', [0,1])
set(h, 'yticklabel', [{'0'}, {'3'}])


elseif strcmp(ParameterSweep, 'NonUniformRing')

x = 0.2 * ([1:22] - 1);
y = 3 * 0.000075 * ([1:22] - 1) ; 

subplot(3,4,9)
Sharpening(StableFP == 0) = nan;
h = imagesc(x, y, Sharpening)
h.AlphaData = ones(size(h.CData)); 
h.AlphaData(isnan(h.CData)) = 0;
set(gca, 'fontsize', 18)
xlabel('E to I Subnetwork Specificity')
ylabel('E to I Subnetwork Weight')
title('Recurrent Sharpening')
set(gca, 'ydir', 'normal')
colormap(gca, colorMap)
h = colorbar
caxis([0,5])

subplot(3,4,10)
SNRFrac(StableFP == 0) = nan;
h = imagesc(x, y, SNRFrac);
h.AlphaData = ones(size(h.CData)); 
h.AlphaData(isnan(h.CData)) = 0;
set(gca, 'fontsize', 18)
xlabel('E to I Subnetwork Specificity')
ylabel('E to I Subnetwork Weight')
title('Norm. Mode Input SNR')
h = colorbar;
caxis([0,1])
set(gca, 'ydir', 'normal')
subplot(3,4,11)
tau(StableFP == 0) = nan;
h = imagesc(x, y, tau);
h.AlphaData = ones(size(h.CData)); 
h.AlphaData(isnan(h.CData)) = 0;
set(gca, 'fontsize', 18)
xlabel('E to I Subnetwork Specificity')
ylabel('E to I Subnetwork Weight')
title('Mode Time Constant')
h = colorbar;
caxis([0,300])
set(gca, 'ydir', 'normal')


subplot(3,4,12)
hold on
for i=1:size(SNRmode,1)
for j=1:size(SNRmode,2)
if StableFP(i,j)
    clr = hsv2rgb([1,1,min(Sharpening(i,j) / 3,1)]);
    scatter(real(sqrt(real(SNRmode{i,j}))), taumode{i,j}, 100, 'marker', '.', 'markeredgecolor', clr, 'linewidth', 5)
end
end
end


set(gca, 'fontsize', 18)
xlabel('Norm. Mode Input SNR')
ylabel('Mode Time Constant')
axis([0,1,0,300])
    
colormap(gca, colorMap)
h = colorbar
set(h,'fontsize', 18)
title(h,'Sharpening')
set(h, 'ytick', [0,1])
set(h, 'yticklabel', [{'0'}, {'3'}])

end



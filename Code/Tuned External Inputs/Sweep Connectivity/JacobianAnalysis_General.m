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

stimvals = 2*pi * [160,200] / 360;
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

for networktype = 1:2

ParameterSweep = 'WeightsvsInputs';

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

elseif strcmp(ParameterSweep, 'EvsIWeights')

kEE = (nEI-1) / stepE;
kIE = (nIE-1) / stepI; 
kEI = 0.4; % I to E concentration
kII = 0.0; 

II_FF_area = IE_FF_area * (0);

elseif strcmp(ParameterSweep, 'WeightsvsInputs')

if networktype == 1    
 
% broad recurrent (sweep inputs)

kEE = 2;
kIE = 0.1; 
kEI = 0.3; % I to E concentration
kII = 0; 

II_FF_area = 0;

kE_FF = 0.1 * nIE;
IE_FF_amp = 0.005 * nEI;
IE_FF_area = IE_FF_amp * pi * besseli(0,kE_FF);

elseif networktype == 2

% narrow recurrent (sweep inputs)    
    
kEE = 3;
kIE = 0.4; 
kEI = 0.4; % I to E concentration
kII = 0; 

II_FF_area = 0;

kE_FF = 0.1 * nIE;
IE_FF_amp = 0.005 * nEI;
IE_FF_area = IE_FF_amp * pi * besseli(0,kE_FF);

end


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
        
        SNRFrac(nEI,nIE) = sqrt(max(real(SNRmode{nEI,nIE})))
        tau(nEI,nIE) = taumode{nEI,nIE}(find(real(SNRmode{nEI,nIE}) == max(real(SNRmode{nEI,nIE})),1))
       
        
    end
end


%% for random search

% find cases where all kXY grow and where SNRFrac increases but tau doesn't
% change much

% plot all networks

% hold on
% for i=1:size(SNRmode,1)
% for j=1:size(SNRmode,2)
% if min(taumode{i,j}) > 0
% scatter(real(sqrt(real(SNRmode{i,j}))), taumode{i,j})
% end
% end
% end


%% Find suitable pairs

SNRFrac0{networktype} = SNRFrac;
tau0{networktype} = tau;

end

tauthresh = 150;
tauratio = tau0{2}(:) ./ tau0{1}(:)';
taumean = (tau0{2}(:) + tau0{1}(:)')/2;
SNRFracratio = SNRFrac0{2}(:) ./ SNRFrac0{1}(:)';

[I,J] = find(and(SNRFracratio > 1, and(and(abs(tauratio - 1) < 0.1, abs(SNRFracratio - 1) > 0.5), taumean > tauthresh)));

[pair1_row, pair1_col] = ind2sub([nEI,nIE],I);
[pair2_row, pair2_col] = ind2sub([nEI,nIE],J);

kE_FF_pair1 = pair1_col * 0.1;
kE_FF_pair2 = pair2_col * 0.1;

IE_FF_amp_pair1 = pair1_row * 0.005;
IE_FF_amp_pair2 = pair2_row * 0.005;

% Find pairs in which connectivity was narrowersub
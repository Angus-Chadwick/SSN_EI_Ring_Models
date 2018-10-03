
clear all

% simulation parameters

Nloop = 10;  % number of simulations for each parameter set
Nvals = 11; % 51
Nt = 5000;  % number of timesteps

% initialise variables

RE_covtot = cell(Nvals);
RE = cell([Nloop, 1]);
RI = RE;
RE_std = RE;
RI_std = RI;
RE_cov = RE;
Rtot_cov = RE;

% fixed input parameters

stimvals = 2*pi * [160,200] / 360;
Nstim = length(stimvals);
theta_aE = 0;
theta_aI = 0;
noise = 0;
kE_FF = 0.5;
kE_TD = 0;
kI_TD = 0;
IE_FF_area = 0.005 * 100;
IE_TD_area = 0.0;
II_TD_area = 0.0;
JEI_mean = 0.04;
JIE_mean = 0.04;

kEI = 0.5;

stepEE = 25 /5 ;
stepIE = 50 / 5;

for nEI = 1:Nvals
    nEI
    for nIE = 1:(5*Nvals)
        
q=1;

% create network        

kEE = (nIE-1) / stepEE; % E to E concentration
kIE = (nEI-1) /stepIE; % I to E concentration

NE = 1000;
network = create_network_varyEE(kEI,kIE,JEI_mean,JIE_mean, 15/NE * besseli(0,1), kEE);

NE = network.cells.NE;
NI = network.cells.NI;

% create inputs

theta_s = stimvals(q);

inputs  = create_inputs(theta_s, theta_aE, theta_aI, noise, kE_FF, IE_FF_area, kE_TD, IE_TD_area, kI_TD, II_TD_area, network);

   
% %% simulate
% 
 NoiseModel = 'Add'; 

[rE, rI]       = SimulateNetwork_mod(network, inputs, Nt, NoiseModel);

R0 = [mean(rE(:,(Nt/2):end),2)', mean(rI(:,(Nt/2):end),2)'];


if ~isnan(sum(R0))
        W = [network.connectivity.JEE, -network.connectivity.JEI; network.connectivity.JIE, -network.connectivity.JII];
        Phip = diag(2 * R0.^(1/2));
        T = diag([network.cells.tauE * ones(NE,1); network.cells.tauI * ones(NI,1)]);
        Jstar = inv(T) * (Phip * W - eye(NE+NI));
        
        DerivativeNorm(nEI,nIE) = norm(inv(T) * ( R0' - max(0, W * R0' + [inputs.IE_FF', inputs.II_FF']').^2));
        
        if DerivativeNorm(nEI,nIE) > 0.001
            
            FixedPointFinder;
            R0 = rmin';
            Phip = diag(2 * R0.^(1/2));
            Jstar = inv(T) * (Phip * W - eye(NE+NI));
        
            DerivativeNorm(nEI,nIE) = norm(inv(T) * ( R0' - max(0, W * R0' + [inputs.IE_FF', inputs.II_FF']').^2));
        
        end
        
        [V,D] = eig(Jstar);
        Evals{nEI,nIE} = diag(D);
        
        inputs.noise = 2;
        
        CovInp = inv(T) * Phip * diag([inputs.noise * mean(inputs.IE_FF) * ones(NE,1); inputs.noise/2 * mean(inputs.IE_FF) * ones(NI,1)]) * Phip * inv(T);
        Inp = (Phip * [inputs.IE_FF .* (- kE_FF * sin(inputs.theta_pE - theta_s))'; zeros(NI,1)]);
        optvec =  Inp ./ diag(CovInp); % works only for diagonal input covariance
        optvec(Inp == 0) = 0;
        optvec = optvec / norm(optvec);
        
        
        Vleft = pinv(V);
        for i=1:size(Vleft,1)
            Vleft(i,:) = Vleft(i,:) ./ norm(  Vleft(i,:));
            angle(i) = -abs(180 /pi * acos(Vleft(i,:) * optvec) - 90) + 90;
        end
        
        
        Angles{nEI,nIE} = angle;
        
end

    end
end


kEE = ([1:nIE]-1) / stepEE;
kIE = ([1:nEI]-1) / stepIE;

for i=1:size(Evals,1)
for j=1:size(Evals,2)
if numel(Evals{i,j}) > 0
l(i,j) = max(real(Evals{i,j}));
lim(i,j) = imag(Evals{i,j}(find(real(Evals{i,j}) == max(real(Evals{i,j})),1))); 
a(i,j) = Angles{i,j}(find(real(Evals{i,j}) == l(i,j),1));
amin(i,j) = min(Angles{i,j});
if ~isnan(amin(i,j))
l_amin(i,j) = real(Evals{i,j}(find(Angles{i,j} == amin(i,j),1))); 
else l_amin(i,j) = nan;
end
else l(i,j) = nan;
end
end
end

subplot(1,3,1)
imagesc(kEE(1:size(l,2)), kIE(1:size(l,1)), (l < 0) + 2 * isnan(l))
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
xlabel('Tuning of E to E weights (k_{EE})')
ylabel('Tuning of E to I weights (k_{IE})')
subplot(1,3,2)
imagesc(kEE(1:size(l,2)), kIE(1:size(l,1)), l)
h = colorbar;
set(gca, 'clim', [-0.1,0.1])
set(gca, 'fontsize', 18)
set(gca, 'ydir', 'normal')
xlabel('Tuning of E to E weights (k_{EE})')
ylabel('Tuning of E to I weights (k_{IE})')
title(h, 'max(real(\lambda))')
subplot(1,3,3)
imagesc(kEE(1:size(l,2)), kIE(1:size(l,1)),a)
h = colorbar;
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
xlabel('Tuning of E to E weights (k_{EE})')
ylabel('Tuning of E to I weights (k_{IE})')
title(h, 'angle')

figure 

subplot(1,3,1)
imagesc(kEE(1:size(l,2)), kIE(1:size(l,1)), (l < 0) + 2 * isnan(l))
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
xlabel('Tuning of E to E weights (k_{EE})')
ylabel('Tuning of E to I weights (k_{IE})')
subplot(1,3,2)
imagesc(kEE(1:size(l,2)), kIE(1:size(l,1)), l_amin)
h = colorbar;
set(gca, 'clim', [-0.1,0.1])
set(gca, 'fontsize', 18)
set(gca, 'ydir', 'normal')
xlabel('Tuning of E to E weights (k_{EE})')
ylabel('Tuning of E to I weights (k_{IE})')
title(h, 'real(\lambda)')
subplot(1,3,3)
imagesc(kEE(1:size(l,2)), kIE(1:size(l,1)),amin)
h = colorbar;
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
xlabel('Tuning of E to E weights (k_{EE})')
ylabel('Tuning of E to I weights (k_{IE})')
title(h, 'min angle')


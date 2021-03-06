 
% Generate tuning curves, standard devation, and selectivity for three
% networks. Output also the external and recurrent currents received by
% each cell. 

clear all

Match_SI = 0; % choose whether to adjust average E<->I weights to match selectivity across network types 

% simulation parameters

Nloop = 10;  % number of simulations for each parameter set
Nt = 10000;  % number of timesteps

% initialise variables

RE = cell([Nloop, 1]);
RI = RE;
RE_std = RE;
RI_std = RI;
RE_cov = RE;
Rtot_cov = RE;

% fixed network parameters 

JEI_mean = 0.04 ;  
JIE_mean = 0.04 ;

% fixed input parameters

theta_aE = 0;
theta_aI = 0;

%% TRY REDUCING NOISE

noise = 0;
%%

kE_FF = 0.5;
kE_TD = 0;
kI_TD = 0;
IE_FF_area = 0.005 * 100;
IE_TD_area = 0.0;
II_TD_area = 0.0;

stimvals = 2*pi * [160, 200] / 360;
Nstim = length(stimvals); 
       
%%

for m=1:2   

    for q=1:Nstim
    
if m==1 % pre

kEE = 2.72;    
kIE = 0.5; 
kEI = 0.5;

elseif m==2 % post

<<<<<<< HEAD
kEE = 3;    
kIE = 0.4; 
=======
kEE = 2.725;
kIE = 0.5;
>>>>>>> f8580b137fe86075998e595c1f43fa95a541b74c
kEI = 0.5;

end
    
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
% 
% 
% parfor n=1:Nloop
% 
%             [rE, rI]       = SimulateNetwork_mod(network, inputs, Nt, NoiseModel);
% 
%     
%     RE{n}(:,q) = mean(rE(:,5000:end),2);
%     RI{n}(:,q) = mean(rI(:,5000:end),2);
%     
%     RE_std{n}(:,q) = std(rE(:,5000:end),[],2);
%     RI_std{n}(:,q) = std(rI(:,5000:end),[],2);
% 
% 
%     RE_cov{n}(:,:,q) = cov([rE(:,5000:end)]')
% 
% end
% 
%    
% 
% 
%     end
% 
% % Analyse results
% 
% RE0 = mean(cat(3,RE{:}),3);
% RI0 = mean(cat(3,RI{:}),3);
% 
% RE0_std = mean(cat(3,RE_std{:}),3);
% RI0_std = mean(cat(3,RI_std{:}),3);
% 
% RE0_cov = mean(cat(4,RE_cov{:}),4);
% 
% 
% mmE = zeros(Nstim);
% mmI = zeros(Nstim);
% m0E = zeros(Nstim);
% m0I = zeros(Nstim);
% 
% iv1 = 1; iv2 = 2;
% 
%         RE0_std_vert{m} = squeeze(RE0_std(:,iv1));
%         RE0_std_ang{m} = squeeze(RE0_std(:,iv2));
% 
%         RE0_vert{m} = squeeze(RE0(:,iv1));
%         RE0_ang{m} = squeeze(RE0(:,iv2));
%      
%         RI0_std_vert{m} = squeeze(RI0_std(:,iv1));
%         RI0_std_ang{m} = squeeze(RI0_std(:,iv2));
% 
%         RI0_vert{m} = squeeze(RI0(:,iv1));
%         RI0_ang{m} = squeeze(RI0(:,iv2));
% 
%         % E selectivity
%         
%         dmu_0E{m} = RE0_ang{m} - RE0_vert{m};
%         
%         poolvar_0E{m} = (0.5* (RE0_std_ang{m}.^2 + RE0_std_vert{m}.^2));
% 
%         SI0E{m} = dmu_0E{m} ./ sqrt(poolvar_0E{m});
%         
%         % I selectivity
%         
%         dmu_0I{m} = RI0_ang{m} - RI0_vert{m};
%         
%         poolvar_0I{m} = (0.5 * (RI0_std_ang{m}.^2 + RI0_std_vert{m}.^2));
% 
%         SI0I{m} = dmu_0I{m} ./ sqrt(poolvar_0I{m});
% 
%         RE_covtot1{m} = RE0_cov(:,:,1);
%         RE_covtot2{m} = RE0_cov(:,:,2);
%        
%         ExternalE_Input{m} = inputs.IE_FF;
%         RecurrentEE_Input{m} = network.connectivity.JEE * RE0_ang{m};
%         RecurrentEI_Input{m} = network.connectivity.JEI * RI0_ang{m};
%     
%         ExternalI_Input{m} = inputs.II_FF;
%         RecurrentIE_Input{m} = network.connectivity.JIE * RE0_ang{m};
%         RecurrentII_Input{m} = network.connectivity.JII * RI0_ang{m};
%         
%         SItot_E_epsridge(m) = squeeze(RE0(:,1) - RE0(:,2))' * inv(0.5 * (RE_covtot1{m} + RE_covtot2{m} + 2 * eps * eye(1000))) * squeeze(RE0(:,1) - RE0(:,2)); 
% 
%         % compute network jacobian 
%         
% 
%         R0       = [RE0_ang{m}', RI0_ang{m}'];
%         
    end
<<<<<<< HEAD

% Analyse results

RE0 = mean(cat(3,RE{:}),3);
RI0 = mean(cat(3,RI{:}),3);

RE0_std = mean(cat(3,RE_std{:}),3);
RI0_std = mean(cat(3,RI_std{:}),3);

RE0_cov = mean(cat(4,RE_cov{:}),4);


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
       
        ExternalE_Input{m} = inputs.IE_FF;
        RecurrentEE_Input{m} = network.connectivity.JEE * RE0_ang{m};
        RecurrentEI_Input{m} = network.connectivity.JEI * RI0_ang{m};
    
        ExternalI_Input{m} = inputs.II_FF;
        RecurrentIE_Input{m} = network.connectivity.JIE * RE0_ang{m};
        RecurrentII_Input{m} = network.connectivity.JII * RI0_ang{m};
        
        SItot_E_epsridge(m) = squeeze(RE0(:,1) - RE0(:,2))' * inv(0.5 * (RE_covtot1{m} + RE_covtot2{m} + 2 * eps * eye(1000))) * squeeze(RE0(:,1) - RE0(:,2)); 

        % compute network jacobian 
        
         inputs_nonoise = inputs;
         inputs_nonoise.noise = 0;
         [rE, rI]       = SimulateNetwork_mod(network, inputs_nonoise, Nt, NoiseModel);
        %[rE, rI]       = SimulateNetwork_mod(network, inputs, Nt, NoiseModel);

        
=======
[rE, rI]       = SimulateNetwork_mod(network, inputs, Nt, NoiseModel);
R0 = [mean(rE(:,5000:end),2)', mean(rI(:,5000:end),2)'];
>>>>>>> f8580b137fe86075998e595c1f43fa95a541b74c
        W = [network.connectivity.JEE, -network.connectivity.JEI; network.connectivity.JIE, -network.connectivity.JII];
        Phip = diag(2 * R0.^(1/2));
        T = diag([network.cells.tauE * ones(NE,1); network.cells.tauI * ones(NI,1)]);
        Jstar{m} = inv(T) * (Phip * W - eye(NE+NI));
%         SE = find(mean(rE(:,5000:end),2) > eps);
%         SE(1) = [];
%         SE(end) = [];
%         SI = find(mean(rI(:,5000:end),2) > eps);
%         SI(1) = [];
%         SI(end) = [];
%         S = zeros([NE + NI,1]);
%         S(SE) = 1;
%         S(NE+SI) = 1;Nt

%         S = logical(S);
        
        Jstar{m} = Jstar{m}; 
        [V{m},D{m}] = eig(Jstar{m});
<<<<<<< HEAD
        inputs1  = create_inputs(theta_s, theta_aE, theta_aI, noise, kE_FF, IE_FF_area, kE_TD, IE_TD_area, kI_TD, II_TD_area, network);
        inputs2  = create_inputs(theta_s * 1.01, theta_aE, theta_aI, noise, kE_FF, IE_FF_area, kE_TD, IE_TD_area, kI_TD, II_TD_area, network);

        u{m} = inv(T) * Phip * [inputs1.IE_FF - inputs2.IE_FF; inputs1.II_FF - inputs2.II_FF];
        Vinv{m} = inv(V{m});
        
        tau{m} = -1./diag(D{m});
        I{m} = find(imag(tau{m}) == 0);
        
        
            for j=1:length(I{m})
                Vinv{m}(I{m}(j),:) = Vinv{m}(I{m}(j),:) ./ norm(Vinv{m}(I{m}(j),:));
            end
            Sigma_eta = diag([ones([1000,1]); 0.5 * ones([200,1])]);
            optvec{m} = inv(inv(T) * Phip) * inv(Sigma_eta) * u{m};
        scatter(tau{m}(I{m}), abs(Vinv{m}(I{m},:) * optvec{m}) ./ norm(optvec{m}));
        hold on
=======
        D{m} = diag(D{m});
        [~,I{m}] = sort(real(D{m}), 'descend');
        D{m} = D{m}(I{m});
        V{m} = V{m}(:,I{m});
        R{m} = find(imag(D{m}) == 0);
% %         hold on
% %         plot(D{m}(R{m}));
%               
        optvec{m} =  inv(Phip) * T * diag(1./[inputs.noise * mean(inputs.IE_FF) * ones(NE,1); inputs.noise/2 * mean(inputs.IE_FF) * ones(NI,1)]) * T * [inputs.IE_FF .* (- kE_FF * sin(inputs.theta_pE - theta_s))'; zeros(NI,1)];
        optvec{m} = optvec{m} / norm(optvec{m});
        
        Vleft{m} = pinv(V{m});
        for i=1:length(R{m})
            Vleft{m}(R{m}(i),:) = Vleft{m}(R{m}(i),:) ./ norm(  Vleft{m}(R{m}(i),:));
            angle{m}(R{m}(i)) = -abs(180 /pi * acos(Vleft{m}(R{m}(i),:) * optvec{m}) - 90) + 90;
        end
%                 
>>>>>>> f8580b137fe86075998e595c1f43fa95a541b74c
end





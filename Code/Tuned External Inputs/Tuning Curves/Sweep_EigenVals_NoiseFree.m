
clear all

% simulation parameters

Nvals = 11; % 51
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

stimvals = 2*pi * [160,200] / 360;
Nstim = length(stimvals);
theta_aE = 0;
theta_aI = 0;
noise = 0;
kE_FF = 0.5; % standard value = 0.5
kE_TD = 0;
kI_TD = 0;
IE_FF_area = 0.005 * 100;
IE_TD_area = 0.0;
II_TD_area = 0.0;
JEI_mean = 0.04;  % standard values
JIE_mean = 0.04;

kEI = 0.4; % standard value 0.5 

stepEE = 25 /5;
stepIE = 50 /5;

for nEI = 1:(2*Nvals)
    nEI
    for nIE = 1:(2*Nvals)
   
q=1;

% create network        
%% 

kEE = (nIE-1) / stepEE; % E to E concentration
kIE = (nEI-1) / stepIE; % I to E concentration

NE = 1000;
network = create_network_varyEE(kEI,kIE,JEI_mean,JIE_mean, 15/NE * besseli(0,1), kEE);

NE = network.cells.NE;
NI = network.cells.NI;

% create inputs

theta_s = stimvals(q);
theta_aI = theta_s;
%kI_TD = kE_FF;
%II_TD_area = 0.1 * IE_FF_area;

inputs  = create_inputs(theta_s, theta_aE, theta_aI, noise, kE_FF, IE_FF_area, kE_TD, IE_TD_area, kI_TD, II_TD_area, network);

   
% %% simulate
 
 NoiseModel = 'Add'; 

[rE, rI]       = SimulateNetwork_mod(network, inputs, Nt, NoiseModel);

StableSim(nEI,nIE) = ~isnan(sum(rE(:)));

if StableSim(nEI,nIE) R0 = [mean(rE(:,(Nt/2):end),2)', mean(rI(:,(Nt/2):end),2)'];else R0 = FixedPoint{nEI,nIE-1};end  % set initial guess
    
            
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
        
            
            
         % NOTE: VERY CAREFULLY CHECK EQUATIONS FOR J, optvec, and Vleft,
         % as T enters differently into different eigenvectors, and have
         % dual basis too...If only E cells count, should be ok, but need
         % to check scaling with time constant and calculation of left dual
         % eigenvector
        
        [Vleft,D] = eig(Jstar');  % try computing left eigenvectors and eigenvalues simply by taking transpose of matrix, no need to invert
        Vleft = Vleft';
        Evals{nEI,nIE} = diag(D);        
        
        inputs.noise = 2;
        Inp = ([inputs.IE_FF .* (- kE_FF * sin(inputs.theta_pE - theta_s))'; zeros(NI,1)]);

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
         
%         Vleft = pinv(V);
%         for i=1:size(Vleft,1)
%             Vleft(i,:) = Vleft(i,:) ./ norm(  Vleft(i,:));
%             angle(i) = -abs(180 /pi * acos(Vleft(i,:) * optvec) - 90) + 90;
%         end

         for i=1:size(Vleft,1)
             Vleft(i,:) = Vleft(i,:);
             angle(i) = -abs(180 /pi * acos(Vleft(i,:)/norm(Vleft(i,1:1000)) * optvec) - 90) + 90;
         end     
%NOTE: taking angle of normed vectors doesn't make sense. The real maths
%tell us we want the dot product of un-normed vectors (Vleft can be
%normed). We can rectify this by simply taking the E part of the e-vec.

        Angles{nEI,nIE} = angle;
        
        
        % compute SNR of best-aligned mode
       % p = find(angle == min(angle),1);
        p = find(max(real(diag(D))),1);
        SNRmode{nEI,nIE} = (Vleft * Inp).^2 ./ diag(Vleft * CovInp * Vleft');
        taumode{nEI,nIE} = -1./real(diag(D));        
        
        
    end
end

kEE = [0:(nIE-1)] / stepEE;
kIE = [0:(nEI-1)] / stepIE;

for i=1:size(Evals,2)
for j=1:size(Evals,1)
M(j,i) = max(real(Evals{j,i}));
end
end

LocalStableFP = (M < 0 & DerivativeNorm < 0.000001);
LocalUnstableFP = (M > 0 & DerivativeNorm < 0.000001);
NoFP = DerivativeNorm > 0.01;
StableSim;

imagesc(kEE,kIE,LocalStableFP + 2 * LocalUnstableFP + 3 * NoFP)

% 
% % 
% [I,J] = sort(real(Evals{nEI,nIE}), 'descend');
% for i=1:3
%         
% subplot(3,2,2*i-1)
% hold on
% plot(real(V(:,J(i))), 'linewidth', 3)
% plot(imag(V(:,J(i))), 'linewidth', 3)
% xlabel('Cell id')
% ylabel('Right Eigenvector Element')
% legend('Real(v)', 'Imag(v)')
% title(strcat('\lambda = ', num2str(Evals{nEI,nIE}(J(i)))))
% 
% subplot(3,2,2*i)
% hold onsub
% plot(real(Vleft(J(i),:)), 'linewidth', 3)
% plot(imag(Vleft(J(i),:)), 'linewidth', 3)
% xlabel('Cell id')
% ylabel('Left Eigenvector Element')
% legend('Real(v)', 'Imag(v)')
% title(strcat('\lambda = ', num2str(Evals{nEI,nIE}(J(i)))))
% 
% end

modechoice = 'SNR';

if strcmp(modechoice, 'timeconstant')

for i=1:size(Evals,1)
for j=1:size(Evals,2)
P = find(abs(imag(Evals{i,j})) < 10e-10);
lambda_max(i,j) = max(real(Evals{i,j}(P)));
ind(i,j) = find(real(Evals{i,j}) == lambda_max(i,j),1);
angle_max(i,j) = Angles{i,j}(ind(i,j));
end
end

elseif strcmp(modechoice, 'angle')

for i=1:size(Evals,1)
for j=1:size(Evals,2)
P = find(abs(imag(Evals{i,j})) < 10e-10);
angle_max(i,j) = min(Angles{i,j}(P));

if ~isnan(angle_max(i,j))
ind(i,j) = find(Angles{i,j} == angle_max(i,j),1);
lambda_max(i,j) = Evals{i,j}(ind(i,j));
else lambda_max(i,j) = nan; end
end
end

elseif strcmp(modechoice, 'SNR')
    
    SNRInp = Inp' * pinv(CovInp) * Inp;
    
    for i=1:size(taumode,1)
        for j=1:size(taumode,2)
                        
           SNR(i,j) = max(real(SNRmode{i,j}));
           angle_max(i,j) = acos(sqrt(SNR(i,j)/SNRInp)) * 180 / pi;
           tau(i,j) = taumode{i,j}(find(real(SNRmode{i,j}) == SNR(i,j),1));
           lambda_max(i,j) = -1./tau(i,j);
            
%             tau(i,j) = max(taumode{i,j});
%             lambda_max(i,j) = -1./tau(i,j); 
%             SNR(i,j) = SNRmode{i,j}(find(taumode{i,j} == tau(i,j),1));
%             angle_max(i,j) = acos(sqrt(SNR(i,j)/SNRInp)) * 180 / pi;
            
        end
    end

end

figure
subplot(2,2,1)
set(gca, 'fontsize', 18)
hold on
for i=[2,6,11]
scatter(kEE(StableSim(i,:)),-1./real(lambda_max(i,StableSim(i,:))), 'linewidth', 3)
xlabel('E to E tuning (k_{EE})')
ylabel('Longest Time Constant (a.u.)')
legend('k_{IE} = 0', 'k_{IE} = 0.5', 'k_{IE} = 1.0')
end
axis([0,6,0,500])
subplot(2,2,2)
set(gca, 'fontsize', 18)
hold on
for i=[2,6,11]
scatter(kEE(StableSim(i,:)),real(angle_max(i,StableSim(i,:))), 'linewidth', 3)
xlabel('E to E tuning (k_{EE})')
ylabel('Angle of Eigenvector from Input Signal')
legend('k_{IE} = 0', 'k_{IE} = 0.5', 'k_{IE} = 1.0')
end
axis([0,6,0,90])
subplot(2,2,3)
set(gca, 'fontsize', 18)
hold on
for i=[6,11,16]
scatter(kIE(StableSim(:,i)),-1./real(lambda_max(StableSim(:,i),i)), 'linewidth', 3)
xlabel('E to I tuning (k_{IE})')
ylabel('Longest Time Constant (a.u.)')
legend('k_{EE} = 1', 'k_{EE} = 2', 'k_{EE} = 3')

end
axis([0,1,0,500])
subplot(2,2,4)
set(gca, 'fontsize', 18)
hold on
for i=[6,11,16]
scatter(kIE(StableSim(:,i)),real(angle_max(StableSim(:,i),i)), 'linewidth', 3)
xlabel('E to I tuning (k_{EI})')
ylabel('Angle of Eigenvector from Input Signal')
legend('k_{EE} = 1', 'k_{EE} = 2', 'k_{EE} = 3')
end
axis([0,1,0,90])



% 
% 
% kEE = ([1:nIE]-1) / stepEE;
% kIE = ([1:nEI]-1) / stepIE;
% 
% for i=1:size(Evals,1)
% for j=1:size(Evals,2)
% if numel(Evals{i,j}) > 0
% l(i,j) = max(real(Evals{i,j}));
% lim(i,j) = imag(Evals{i,j}(find(real(Evals{i,j}) == max(real(Evals{i,j})),1))); 
% a(i,j) = Angles{i,j}(find(real(Evals{i,j}) == l(i,j),1));
% amin(i,j) = min(Angles{i,j});
% if ~isnan(amin(i,j))
% l_amin(i,j) = real(Evals{i,j}(find(Angles{i,j} == amin(i,j),1))); 
% else l_amin(i,j) = nan;
% end
% else l(i,j) = nan;
% end
% end
% end
% 
% subplot(1,3,1)
% imagesc(kEE(1:size(l,2)), kIE(1:size(l,1)), (l < 0) + 2 * isnan(l))
% set(gca, 'ydir', 'normal')
% set(gca, 'fontsize', 18)
% xlabel('Tuning of E to E weights (k_{EE})')
% ylabel('Tuning of E to I weights (k_{IE})')
% subplot(1,3,2)
% imagesc(kEE(1:size(l,2)), kIE(1:size(l,1)), l)
% h = colorbar;
% set(gca, 'clim', [-0.1,0.1])
% set(gca, 'fontsize', 18)
% set(gca, 'ydir', 'normal')
% xlabel('Tuning of E to E weights (k_{EE})')
% ylabel('Tuning of E to I weights (k_{IE})')
% title(h, 'max(real(\lambda))')
% subplot(1,3,3)
% imagesc(kEE(1:size(l,2)), kIE(1:size(l,1)),a)
% h = colorbar;
% set(gca, 'ydir', 'normal')
% set(gca, 'fontsize', 18)
% xlabel('Tuning of E to E weights (k_{EE})')
% ylabel('Tuning of E to I weights (k_{IE})')
% title(h, 'angle')
% 
% figure 
% 
% subplot(1,3,1)
% imagesc(kEE(1:size(l,2)), kIE(1:size(l,1)), (l < 0) + 2 * isnan(l))
% set(gca, 'ydir', 'normal')
% set(gca, 'fontsize', 18)
% xlabel('Tuning of E to E weights (k_{EE})')
% ylabel('Tuning of E to I weights (k_{IE})')
% subplot(1,3,2)
% imagesc(kEE(1:size(l,2)), kIE(1:size(l,1)), l_amin)
% h = colorbar;
% set(gca, 'clim', [-0.1,0.1])
% set(gca, 'fontsize', 18)
% set(gca, 'ydir', 'normal')
% xlabel('Tuning of E to E weights (k_{EE})')
% ylabel('Tuning of E to I weights (k_{IE})')
% title(h, 'real(\lambda)')
% subplot(1,3,3)
% imagesc(kEE(1:size(l,2)), kIE(1:size(l,1)),amin)
% h = colorbar;
% set(gca, 'ydir', 'normal')
% set(gca, 'fontsize', 18)
% xlabel('Tuning of E to E weights (k_{EE})')
% ylabel('Tuning of E to I weights (k_{IE})')
% title(h, 'min angle')


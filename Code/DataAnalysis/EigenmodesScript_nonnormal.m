
if task='learning'

elseif task='switching'
    
xpre = xirrel;
xpost = xrel;
residualpre = residualirrel;
residualpost = residualrel;

end
    
cond = 'differential'

for rix=1:length(RXS)

    % compute signal vectors
    
          MeanVpre = mean(xpre{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9):(length(CELLLAB_ALL{RXS(rix)}) + 17)), 2);
          MeanApre = mean(xpre{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9 + 17):(length(CELLLAB_ALL{RXS(rix)}) + 17 + 17)), 2);

          MeanVpost = mean(xpost{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9):(length(CELLLAB_ALL{RXS(rix)}) + 17)), 2);
          MeanApost = mean(xpost{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9 + 17):(length(CELLLAB_ALL{RXS(rix)}) + 17 + 17)), 2);
  
    
          if strmatch(cond, 'differential')

            dIn_pre{rix} = MeanVpre - MeanApre;
            dIn_post{rix} = MeanVpost - MeanApost;
    
          elseif strmatch(cond, 'Vert') 
              
            dIn_pre{rix} = MeanVpre;
            dIn_post{rix} = MeanVpost;
            
          elseif strmatch(cond, 'Ang') 
              
            dIn_pre{rix} = MeanApre;
            dIn_post{rix} = MeanApost;
            
         
          elseif strmatch(cond, 'Shuffled')

            dIn_pre{rix} = MeanVpre - MeanApre;
            dIn_post{rix} = MeanVpost - MeanApost;
            
            dIn_pre{rix} = dIn_pre{rix}(randperm(length(dIn_pre{rix})));
            dIn_post{rix} = dIn_post{rix}(randperm(length(dIn_post{rix})));
            
          end
    
            
            Covres_pre{rix} = cov(residualpre{RXS(rix)}');
            Covres_post{rix} = cov(residualpost{RXS(rix)}');
            Covresinv_pre{rix} = inv(Covres_pre{rix});
            Covresinv_post{rix} = inv(Covres_post{rix});

            InputInformation_pre(rix) = dIn_pre{rix}' * Covresinv_pre{rix} * dIn_pre{rix};
            InputInformation_post(rix) = dIn_post{rix}' * Covresinv_post{rix} * dIn_post{rix};

                        
    % compute eigendecomposition of weight matrix
          
    Apre{rix} = xpre{RXS(rix)}(:, 1:size(xpre{RXS(rix)},1));
    Apost{rix} = xpost{RXS(rix)}(:, 1:size(xpost{RXS(rix)},1));
    
    [Vpre{rix}, Dpre{rix}] = eig(Apre{rix});
    [Vpost{rix}, Dpost{rix}] = eig(Apost{rix});
        
    Vinv_pre{rix} = inv(Vpre{rix});
    Vinv_post{rix} = inv(Vpost{rix});
    
    
    SNRfrac_pre{rix}  = abs(Vinv_pre{rix}  * dIn_pre{rix})  ./ sqrt(diag(Vinv_pre{rix}  * Covres_pre{rix}  * Vinv_pre{rix}'))  / sqrt(InputInformation_pre(rix));
    SNRfrac_post{rix} = abs(Vinv_post{rix} * dIn_post{rix}) ./ sqrt(diag(Vinv_post{rix} * Covres_post{rix} * Vinv_post{rix}')) / sqrt(InputInformation_post(rix));

    % compute non-normality metrics
    
%     DModes_pre{rix} = Vinv_pre{rix} * Covres_pre{rix} * Vinv_pre{rix}';
%     DModes_post{rix} = Vinv_post{rix} * Covres_post{rix} * Vinv_post{rix}';
    
%     DcorrModes_pre{rix}  = DModes_pre{rix}./sqrt(diag(DModes_pre{rix}) * diag(DModes_pre{rix})') + tril(nan(size(DModes_pre{rix})),0);
%     DcorrModes_post{rix} = DModes_post{rix}./sqrt(diag(DModes_post{rix}) * diag(DModes_post{rix})') + tril(nan(size(DModes_post{rix})),0);
  
    AlignmentModes_pre{rix} = acos(Vpre{rix}' * Vpre{rix});
    AlignmentModes_post{rix} = acos(Vpost{rix}' * Vpost{rix});

    % subselect real eigenmodes
    
    I = find(imag(diag(Dpre{rix})) == 0);
    D = diag(Dpre{rix});
    RealEig_pre{rix} = real(D(I));
    SNRfrac_pre{rix} = real(SNRfrac_pre{rix}(I));
 %   DcorrModes_pre{rix} = real(DcorrModes_pre{rix}(I,I));
    AlignmentModes_pre{rix} = real(AlignmentModes_pre{rix}(I,I)) + tril(nan(length(RealEig_pre{rix})),0);
    
    
    I = find(imag(diag(Dpost{rix})) == 0);
    D = diag(Dpost{rix});
    RealEig_post{rix} = real(D(I));
    SNRfrac_post{rix} = real(SNRfrac_post{rix}(I));
%    DcorrModes_post{rix} = real(DcorrModes_post{rix}(I,I));
    AlignmentModes_post{rix} = real(AlignmentModes_post{rix}(I,I)) + tril(nan(length(RealEig_post{rix})),0);
    
%     RealEig_Fracs_pre{rix}  = bsxfun(@times, RealEig_pre{rix}, 1./RealEig_pre{rix}') + tril(nan(length(RealEig_pre{rix})),0);
%     RealEig_Fracs_post{rix} = bsxfun(@times, RealEig_post{rix}, 1./RealEig_post{rix}') + tril(nan(length(RealEig_post{rix})),0);      


    RealEig_Diffs_pre{rix}  = bsxfun(@minus, RealEig_pre{rix}, RealEig_pre{rix}') + tril(nan(length(RealEig_pre{rix})),0);
    RealEig_Diffs_post{rix} = bsxfun(@minus, RealEig_post{rix}, RealEig_post{rix}') + tril(nan(length(RealEig_post{rix})),0);      


    

%     DcorrModes_pre{rix} = DcorrModes_pre{rix}(:);
%     DcorrModes_post{rix} = DcorrModes_post{rix}(:);
%     RealEig_Fracs_pre{rix} = RealEig_Fracs_pre{rix}(:);
%     RealEig_Fracs_post{rix} = RealEig_Fracs_post{rix}(:);
% 
%     AlignmentModes_pre{rix} = AlignmentModes_pre{rix}(:);
%     AlignmentModes_post{rix} = AlignmentModes_post{rix}(:);
%     RealEig_Diffs_pre{rix} = RealEig_Diffs_pre{rix}(:);
%     RealEig_Diffs_post{rix} = RealEig_Diffs_post{rix}(:);


%     NonNormalAmp_pre{rix} = 1 - 4 * RealEig_Fracs_pre{rix} ./ (1+RealEig_Fracs_pre{rix}).^2 .* DcorrModes_pre{rix}.^2;
%     NonNormalAmp_post{rix} = 1 - 4 * RealEig_Fracs_post{rix} ./ (1+RealEig_Fracs_post{rix}).^2 .* DcorrModes_post{rix}.^2;
    
%      NonNormalAmp_pre{rix} = RealEig_Diffs_pre{rix} .* cos(AlignmentModes_pre{rix}) ./ sin(AlignmentModes_pre{rix}) ;
%      NonNormalAmp_post{rix} = RealEig_Diffs_post{rix} .* cos(AlignmentModes_post{rix}) ./ sin(AlignmentModes_post{rix}) ;

     NonNormalAmp_pre{rix} = -abs(RealEig_Diffs_pre{rix}) .* cos(AlignmentModes_pre{rix})  ;
     NonNormalAmp_post{rix} = -abs(RealEig_Diffs_post{rix}) .* cos(AlignmentModes_post{rix})  ;
     


end


%% Try to plot measures of multiple defectiveness

for rix=1:length(RXS)
subplot(2,1,1)
hold on
X = -abs(RealEig_Diffs_pre{rix}); % 
Y = abs(cos(AlignmentModes_pre{rix}));
Xtot_pre{rix} = nansum(1./X) + nansum(1./X');
Ytot_pre{rix} = nansum(Y) + nansum(Y');
scatter(Xtot_pre{rix} , (Ytot_pre{rix} ))  % X is similarity of eigenvalues


subplot(2,1,2)
hold on
X = -abs(RealEig_Diffs_post{rix});
Y = abs(cos(AlignmentModes_post{rix}));
Xtot_post{rix} = nansum(1./X) + nansum(1./X');
Ytot_post{rix} = nansum(Y) + nansum(Y');
scatter(Xtot_post{rix}, (Ytot_post{rix}) )
end


%%
    

plot_inputinf = 0;
if plot_inputinf
scatter(InputInformation_pre, InputInformation_post, 100, 'linewidth', 3)
hold on
plot([0,100], [0,100], 'linewidth', 3)
set(gca, 'fontsize', 18)
xlabel('Input Information (pre)')
ylabel('Input Information (post)')
title('Change in Input Information with Learning')
p=signrank(InputInformation_pre, InputInformation_post);
legend(strcat('p=', num2str(p)))
box on
end

rix = 2;
    
% combine animals
% 
 RealEig_pre_tot = vertcat(RealEig_pre{:});
% SNRfrac_pre_tot = vertcat(SNRfrac_pre{:});
 RealEig_post_tot = vertcat(RealEig_post{:});
% SNRfrac_post_tot = vertcat(SNRfrac_post{:});
% 
% RealEig_Fracs_pre_tot =  vertcat(RealEig_Fracs_pre{:});
% RealEig_Fracs_post_tot = vertcat(RealEig_Fracs_post{:});
% DcorrModes_pre_tot = vertcat(DcorrModes_pre{:});
% DcorrModes_post_tot = vertcat(DcorrModes_post{:});
%  
% NonNormalAmp_pre_tot = 1 - 4 * RealEig_Fracs_pre_tot ./ (1+RealEig_Fracs_pre_tot).^2 .* DcorrModes_pre_tot.^2;
% NonNormalAmp_post_tot = 1 - 4 * RealEig_Fracs_post_tot ./ (1+RealEig_Fracs_post_tot).^2 .* DcorrModes_post_tot.^2;
% 
% 
% NonNormalAmp_pre_tot(isnan(NonNormalAmp_pre_tot)) = [];
% NonNormalAmp_post_tot(isnan(NonNormalAmp_post_tot)) = [];

RealEig_Diffs_pre_tot = vertcat(RealEig_Diffs_pre{:});
RealEig_Diffs_post_tot = vertcat(RealEig_Diffs_post{:});
NonNormalAmp_pre_tot = vertcat(NonNormalAmp_pre{:});
NonNormalAmp_post_tot = vertcat(NonNormalAmp_post{:});
% select eigenmodes within physical bounds (stable and not excessively long time constant)

Ipre = find(and(0 < RealEig_pre_tot + 1, RealEig_pre_tot + 1 < 0.99));
Ipost = find(and(0 < RealEig_post_tot + 1, RealEig_post_tot + 1 < 0.99));

% convert discrete eigenvalues to continuous time constants

taupre_tot = -125./log(RealEig_pre_tot(Ipre)+1);
taupost_tot = -125./log(RealEig_post_tot(Ipost)+1);

% take moving average windows of eigenmodes (time windows and angle windows)

binwidth = 100;
binstep = 0.25 * binwidth;
taubins = binwidth:binstep:1400;

SNRfracbin_width = 0.025;
SNRfracbin_step = 0.25 * SNRfracbin_width;
SNRbins = SNRfracbin_width:SNRfracbin_step:0.25; 

clear SNRfracbin_pre SNRfracbin_post SNRfracSEMbin_pre SNRfracSEMbin_post p1 p2

for i=1:length(taubins)
      
    eigs_pre = and(taupre_tot >= taubins(i) - binwidth, taupre_tot < taubins(i) + binwidth);
    eigs_post = and(taupost_tot >= taubins(i) - binwidth, taupost_tot < taubins(i) + binwidth);
    Neigs_pre = sum(eigs_pre);
    Neigs_post = sum(eigs_post);  
   
    % compute mean and variance of SNRfracs
   
 
    SNRfracbin_pre(i) = mean(SNRfrac_pre_tot(Ipre( eigs_pre )));
    SNRfracbin_post(i) = mean(SNRfrac_post_tot(Ipost( eigs_post )));

    SNRfracSEMbin_pre(i) = std(SNRfrac_pre_tot(Ipre( eigs_pre ))) ./ sqrt(Neigs_pre);
    SNRfracSEMbin_post(i) = std(SNRfrac_post_tot(Ipost( eigs_post )))  ./ sqrt(Neigs_post);
    
    
end

    clear taubin_pre taubin_post taubin_pre_sem taubin_post_sem
    
    for i=1:length(SNRbins)
        
       eigs_pre = and( SNRfrac_pre_tot(Ipre) >= SNRbins(i) - SNRfracbin_width, SNRfrac_pre_tot(Ipre) < SNRbins(i) + SNRfracbin_width);
       eigs_post = and( SNRfrac_post_tot(Ipost) >= SNRbins(i) - SNRfracbin_width,  SNRfrac_post_tot(Ipost) < SNRbins(i) + SNRfracbin_width);
       Neigs_pre = sum(eigs_pre);
       Neigs_post = sum(eigs_post);       
       
       
    taubin_pre(i) = mean(taupre_tot(eigs_pre));
    taubin_post(i) = mean(taupost_tot(eigs_post));

    taubin_pre_sem(i) = std(taupre_tot( eigs_pre)) / sqrt(Neigs_pre);
    taubin_post_sem(i) = std(taupost_tot( eigs_post)) / sqrt(Neigs_post);
    
  %  p2(i) = ranksum(taupre_tot(eigs_pre), taupost_tot(eigs_post));
    
    end
    
   plotrawdata=0;
 
subplot(1,2,1)
shadedErrorBar(taubins, SNRfracbin_pre , SNRfracSEMbin_pre ,'lineprops','k');
hold on
shadedErrorBar(taubins, SNRfracbin_post , SNRfracSEMbin_post,'lineprops','b');
%axis([0,1250,75.5,90])
set(gca, 'fontsize', 18)
xlabel('Time Constant (ms)')
ylabel('Fraction of Feedforward Information')

if plotrawdata
hold on
    scatter(taupre_tot,  SNRfrac_pre_tot(Ipre), 100, 'x', 'markeredgecolor', 'k', 'linewidth', 3)
    scatter(taupost_tot,  SNRfrac_post_tot(Ipost), 100, 'o', 'markeredgecolor', 'b', 'linewidth', 3)       
    
    legend('Pre-Learning' ,'Post-Learning')
    box on
    xlabel('Time constant (ms)')
    ylabel('Fraction of Feedforward Information')
    title('Dynamical Modes')
    set(gca, 'fontsize', 24)
    
end


subplot(1,2,2)
    
shadedErrorBar(SNRbins, taubin_pre, taubin_pre_sem,'lineprops','k');
hold on
shadedErrorBar(SNRbins, taubin_post, taubin_post_sem,'lineprops','b');
set(gca, 'fontsize', 18)
ylabel('Time Constant (ms)')
xlabel('Fraction of Feedforward Information')
%axis([75.5,90,0,1250])


%% smoothed 2d hists

clear p_pre_smx p_post_smx p_pre_smxy p_post_smxy



edges{1} = 0:10:1750;
edges{2} = 0:0.001:0.25;

sigma_t = 100;
sigma_SNRfrac = 0.025;

s1=(exp(-(bsxfun(@minus, edges{1}, taupre_tot)).^2/(2 * sigma_t^2)));
s2=(exp(-(bsxfun(@minus, edges{2}, SNRfrac_pre_tot(Ipre)).^2)/(2*sigma_SNRfrac^2)));
s3=s1' * s2;
p_pre_smxy = s3;
s1=(exp(-(bsxfun(@minus, edges{1}, taupost_tot)).^2/(2 * sigma_t^2)));
s2=(exp(-(bsxfun(@minus, edges{2}, SNRfrac_post_tot(Ipost)).^2)/(2*sigma_SNRfrac^2)));
s3=s1' * s2;
p_post_smxy = s3;

normtype = 'shuffstd'

if strcmp(normtype, 'unnormalised')


figure 
cmax = 0.05;
subplot(3,1,1)
imagesc(edges{1}, edges{2}, p_pre_smxy')
h = colorbar;
%caxis([0,cmax])
%axis([0,1250,76,90])
subplot(3,1,2)
imagesc(edges{1}, edges{2} ,p_post_smxy')
h = colorbar;
%caxis([0,cmax])
%axis([0,1250,76,90])
subplot(3,1,3)
imagesc(edges{1}, edges{2} ,(p_post_smxy' - p_pre_smxy'))
h = colorbar;
%caxis([-cmax,cmax])
%axis([0,1250,76,90])

elseif strcmp(normtype, 'shuffstd')
  
    prc = 2.5;
    
for i=1:length(edges{1})
for j=1:length(edges{2})
normvar_diff(i,j) = std(p_post_smxy_shuff(i,j,:) - p_pre_smxy_shuff(i,j,:));
normvar_pre(i,j) = std(p_pre_smxy_shuff(i,j,:));
normvar_post(i,j) = std(p_post_smxy_shuff(i,j,:));
p_top_prc(i,j) = prctile(p_post_smxy_shuff(i,j,:) - p_pre_smxy_shuff(i,j,:),100 - prc);
p_bottom_prc(i,j) = prctile(p_post_smxy_shuff(i,j,:) - p_pre_smxy_shuff(i,j,:),prc);
end
end



%cmax1 = 10;
%cmax2 = 3;

figure 
subplot(3,1,1)
imagesc(edges{1}, edges{2} , p_pre_smxy' ./ normvar_pre')
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
xlabel('Time constant (ms)')
ylabel('Angle from optimal input SNR (deg)')
h = colorbar;
title(h, 'Normalized count')
%caxis([0,cmax1])
%axis([0,1250,76,90])
subplot(3,1,2)
imagesc(edges{1}, edges{2} ,p_post_smxy' ./ normvar_post')
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
xlabel('Time constant (ms)')
ylabel('Angle from optimal input SNR (deg)')
h = colorbar;
title(h, 'Normalized count')
%caxis([0,cmax1])
%axis([0,1250,76,90])
subplot(3,1,3)
imagesc(edges{1}, edges{2} ,(p_post_smxy' - p_pre_smxy') ./ normvar_diff')
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
xlabel('Time constant (ms)')
ylabel('Angle from optimal input SNR (deg)')
h = colorbar;
title(h, 'Normalized count')
%caxis([-cmax2,cmax2])
%axis([0,1250,76,90])
hold on
X1 = (p_post_smxy' - p_pre_smxy' > p_top_prc');
X2 = (p_post_smxy' - p_pre_smxy' < p_bottom_prc');
contour(edges{1},edges{2},X1, [0,1], 'color', [0.5,0.5,0.5], 'linewidth', 4, 'linestyle', '--');
contour(edges{1},edges{2},X2, [0,1], 'k', 'linewidth', 4, 'linestyle', '--');

end

%% smoothed hists nonnormal

clear p_pre_smx p_post_smx p_pre_smxy p_post_smxy

edges{2} = 0:0.01:1.4;
edges{1} = 0:0.01:0.35;

sigma_1 = 0.05;
sigma_2 = 0.05;

s1=(exp(-(bsxfun(@minus, edges{1}, abs(NonNormalAmp_pre_tot(~isnan(RealEig_Diffs_pre_tot)))) ).^2/(2 * sigma_1^2)));
s2=(exp(-(bsxfun(@minus, edges{2}, abs(RealEig_Diffs_pre_tot(~isnan(RealEig_Diffs_pre_tot)))) .^2)/(2*sigma_2^2)));
s3=s1' * s2;
p_pre_smxy = s3;
s1=(exp(-(bsxfun(@minus, edges{1}, abs(NonNormalAmp_post_tot(~isnan(RealEig_Diffs_post_tot))))) .^2/(2 * sigma_1^2)));
s2=(exp(-(bsxfun(@minus, edges{2}, abs(RealEig_Diffs_post_tot(~isnan(RealEig_Diffs_post_tot)))) .^2)/(2*sigma_2^2)));
s3=s1' * s2;
p_post_smxy = s3;


p_pre_smxy = p_pre_smxy / sum(p_pre_smxy(:));
p_post_smxy = p_post_smxy / sum(p_post_smxy(:));

figure 
cmax = 0.05;
subplot(3,1,1)
imagesc(edges{1}, edges{2}, p_pre_smxy')
h = colorbar;
%caxis([0,cmax])
%axis([0,1250,76,90])
subplot(3,1,2)
imagesc(edges{1}, edges{2} ,p_post_smxy')
h = colorbar;
%caxis([0,cmax])
%axis([0,1250,76,90])
subplot(3,1,3)
imagesc(edges{1}, edges{2} ,(p_post_smxy' - p_pre_smxy'))
h = colorbar;
%caxis([-cmax,cmax])
%axis([0,1250,76,90])

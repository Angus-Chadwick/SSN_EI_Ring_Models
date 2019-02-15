
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
            
            sub
          elseif strmatch(cond, 'Shuffled')

            dIn_pre{rix} = MeanVpre - MeanApre;
            dIn_post{rix} = MeanVpost - MeanApost;
            
            dIn_pre{rix} = dIn_pre{rix}(randperm(length(dIn_pre{rix})));
            dIn_post{rix} = dIn_post{rix}(randperm(length(dIn_post{rix})));
            
          end
    
          % compute SNR vector and information

          SNRvec = 1;
          if SNRvec
          
            Covres_pre{rix} = cov(residualpre{RXS(rix)}');
            Covres_post{rix} = cov(residualpost{RXS(rix)}');
            Covresinv_pre{rix} = inv(Covres_pre{rix});
            Covresinv_post{rix} = inv(Covres_post{rix});

            InputInformation_pre(rix) = dIn_pre{rix}' * Covresinv_pre{rix} * dIn_pre{rix};
            InputInformation_post(rix) = dIn_post{rix}' * Covresinv_post{rix} * dIn_post{rix};

            dIn_pre{rix} = Covresinv_pre{rix} * dIn_pre{rix};
            dIn_post{rix} = Covresinv_post{rix} * dIn_post{rix};
            
          end
    
          % compute eigendecomposition
          
    Apre{rix} = xpre{RXS(rix)}(:, 1:size(xpre{RXS(rix)},1));
    Apost{rix} = xpost{RXS(rix)}(:, 1:size(xpost{RXS(rix)},1));
    
    [Vpre{rix}, Dpre{rix}] = eig(Apre{rix});
    [Vpost{rix}, Dpost{rix}] = eig(Apost{rix});
        
    Vinv_pre{rix} = inv(Vpre{rix});
    Vinv_post{rix} = inv(Vpost{rix});
    
    for i=1:size(Vinv_post{rix},1)
    
        Vinv_pre{rix}(i,:)  = Vinv_pre{rix}(i,:) ./ norm(Vinv_pre{rix}(i,:));
        Vinv_post{rix}(i,:) = Vinv_post{rix}(i,:) ./ norm(Vinv_post{rix}(i,:));
        
    end
   
    % subselect real eigenmodes
    
    I = find(imag(diag(Dpre{rix})) == 0);
    Xpre{rix} = abs(Vinv_pre{rix}(I,:) * dIn_pre{rix}) / norm(dIn_pre{rix});
    D = diag(Dpre{rix});
    Ypre{rix} = D(I);
    
    I = find(imag(diag(Dpost{rix})) == 0);
    Xpost{rix} = abs(Vinv_post{rix}(I,:) * dIn_post{rix}) / norm(dIn_post{rix});
    D = diag(Dpost{rix});
    Ypost{rix} = D(I);
  
   % plot significant eigenvectors vs input SNR
%     
%     l = and(and(0 < Ypost{rix} + 1, Ypost{rix} + 1 < 0.99),and(and(-125./log(Ypost{rix}+1) > 800, -125./log(Ypost{rix}+1) < 1000), and(acos(Xpost{rix}) * 180 / pi > 82, acos(Xpost{rix}) * 180 / pi < 86)));
%   
%     
%     figure
%     if sum(l) > 0
%         p=0;
%         s=find(l);
%     for t=1:length(s)
%         q = (s(t));
%         p=p+1;
%         subplot(1,sum(l),p)
%             hold on
% 
% 
%             
%             [c1{rix}(t), pval{rix}(t)] = corr(dIn_post{rix}(CELLLAB_ALL{RXS(rix)}==1),real(Vinv_pre{rix}(l(q),CELLLAB_ALL{RXS(rix)}==1))');
%             if c1{rix}(t) > 0
%                 
%             c3{rix}(t) = corr(dIn_post{rix}(CELLLAB_ALL{RXS(rix)}==3),real(Vinv_pre{rix}(l(q),CELLLAB_ALL{RXS(rix)}==3))');
%             else
%                 c1{rix}(t) = -c1{rix}(t);
%                 c3{rix}(t) = -corr(dIn_post{rix}(CELLLAB_ALL{RXS(rix)}==3),real(Vinv_pre{rix}(l(q),CELLLAB_ALL{RXS(rix)}==3))');
%             end
%            
%             scatter(dIn_post{rix}(CELLLAB_ALL{RXS(rix)}==1), real(Vinv_pre{rix}(l(q),CELLLAB_ALL{RXS(rix)}==1)), 'r')
%             scatter(dIn_post{rix}(CELLLAB_ALL{RXS(rix)}==3), real(Vinv_pre{rix}(l(q),CELLLAB_ALL{RXS(rix)}==3)), 'g')
%             scatter(dIn_post{rix}(CELLLAB_ALL{RXS(rix)}==4), real(Vinv_pre{rix}(l(q),CELLLAB_ALL{RXS(rix)}==4)), 'b')
%             scatter(dIn_post{rix}(CELLLAB_ALL{RXS(rix)}==5), real(Vinv_pre{rix}(l(q),CELLLAB_ALL{RXS(rix)}==5)), 'k')
% 
%     end
%     end
%     
 end


plot_inputinf = 1;
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
combine_animals = 1;

rix = 2

if combine_animals
    
% combine animals

Xpre_tot = vertcat(Xpre{:});
Ypre_tot = vertcat(Ypre{:});
Xpost_tot = vertcat(Xpost{:});
Ypost_tot = vertcat(Ypost{:});

else 
    
 
Xpre_tot = (Xpre{rix});
Ypre_tot = (Ypre{rix});
Xpost_tot = (Xpost{rix});
Ypost_tot = (Ypost{rix});   


end

% compute angles

thetapre_tot = acos(Xpre_tot);
thetapost_tot = acos(Xpost_tot);


% select eigenmodes within physical bounds (stable and not excessively long time constant)

Ipre = find(and(0 < Ypre_tot + 1, Ypre_tot + 1 < 0.99));
Ipost = find(and(0 < Ypost_tot + 1, Ypost_tot + 1 < 0.99));

% convert discrete eigenvalues to continuous time constants

taupre_tot = -125./log(Ypre_tot(Ipre)+1);
taupost_tot = -125./log(Ypost_tot(Ipost)+1);

% take moving average windows of eigenmodes (time windows and angle windows)

binwidth = 100;
binstep = 0.25 * binwidth;
taubins = binwidth:binstep:1400;

thetabin_width = 1.0;
thetabinstep = 0.25 * thetabin_width;
thetabins = 78:thetabinstep:(90 - thetabin_width); 

clear thetabin_pre thetabin_post thetaSEMbin_pre thetaSEMbin_post p1 p2

for i=1:length(taubins)
      
    eigs_pre = and(taupre_tot >= taubins(i) - binwidth, taupre_tot < taubins(i) + binwidth);
    eigs_post = and(taupost_tot >= taubins(i) - binwidth, taupost_tot < taubins(i) + binwidth);
    Neigs_pre = sum(eigs_pre);
    Neigs_post = sum(eigs_post);  
   
    % compute circular mean and variance of angles
   
    resultantbin_pre(i) = mean(exp(1i * thetapre_tot( Ipre( eigs_pre ))));
    resultantbin_post(i) = mean(exp(1i * thetapost_tot( Ipost( eigs_post ))));

    thetabin_pre(i) = angle(resultantbin_pre(i));
    thetabin_post(i) = angle(resultantbin_post(i));

    thetaSEMbin_pre(i) = sqrt((1 - abs(resultantbin_pre(i)))) ./ sqrt(Neigs_pre);
    thetaSEMbin_post(i) = sqrt((1 - abs(resultantbin_post(i))))  ./ sqrt(sum( and(taupost_tot >= taubins(i) - binwidth, taupost_tot < taubins(i) + binwidth)));
    
    
 %   p1(i) = ranksum(thetapre_tot(Ipre(eigs_pre)), thetapost_tot(Ipost(eigs_post)));

end

    clear taubin_pre taubin_post taubin_pre_sem taubin_post_sem
    
    for i=1:length(thetabins)
        
       eigs_pre = and(180 / pi * thetapre_tot(Ipre) >= thetabins(i) - thetabin_width, 180 / pi * thetapre_tot(Ipre) < thetabins(i) + thetabin_width);
       eigs_post = and(180 / pi * thetapost_tot(Ipost) >= thetabins(i) - thetabin_width, 180 / pi * thetapost_tot(Ipost) < thetabins(i) + thetabin_width);
       Neigs_pre = sum(eigs_pre);
       Neigs_post = sum(eigs_post);       
       
       
    taubin_pre(i) = mean(taupre_tot(eigs_pre));
    taubin_post(i) = mean(taupost_tot(eigs_post));

    taubin_pre_sem(i) = std(taupre_tot( eigs_pre)) / sqrt(Neigs_pre);
    taubin_post_sem(i) = std(taupost_tot( eigs_post)) / sqrt(Neigs_post);
    
  %  p2(i) = ranksum(taupre_tot(eigs_pre), taupost_tot(eigs_post));
    
    end
    
   plotrawdata=1;
 
 subplot(1,2,1)
shadedErrorBar(taubins, thetabin_pre * 180/pi, thetaSEMbin_pre * 180/pi,'lineprops','k');
hold on
shadedErrorBar(taubins, thetabin_post * 180/pi, thetaSEMbin_post * 180/pi,'lineprops','b');
axis([0,1250,75.5,90])
set(gca, 'fontsize', 18)
xlabel('Time Constant (ms)')
ylabel('Angle from Optimal Input SNR')

if plotrawdata
hold on
    scatter(taupre_tot, 180 / pi * thetapre_tot(Ipre), 100, 'x', 'markeredgecolor', 'k', 'linewidth', 3)
    scatter(taupost_tot, 180 / pi * thetapost_tot(Ipost), 100, 'o', 'markeredgecolor', 'b', 'linewidth', 3)       
    
    legend('Pre-Learning' ,'Post-Learning')
    box on
    xlabel('Time constant (ms)')
    ylabel('Angle from optimal input SNR (deg)')
    title('Dynamical Modes')
    set(gca, 'fontsize', 24)
    
end


subplot(1,2,2)
    
shadedErrorBar(thetabins, taubin_pre, taubin_pre_sem,'lineprops','k');
hold on
shadedErrorBar(thetabins, taubin_post, taubin_post_sem,'lineprops','b');
set(gca, 'fontsize', 18)
ylabel('Time Constant (ms)')
xlabel('Angle from Optimal Input SNR')
axis([75.5,90,0,1250])


%% smoothed 2d hists

clear p_pre_smx p_post_smx p_pre_smxy p_post_smxy



edges{1} = 0:10:1750;
edges{2} = (70:0.1:90) * pi / 180;

sigma_t = 100;
sigma_theta = 1 * pi / 180;

s1=(exp(-(bsxfun(@minus, edges{1}, taupre_tot)).^2/(2 * sigma_t^2)));
s2=(exp(-(bsxfun(@minus, edges{2}, thetapre_tot(Ipre)).^2)/(2*sigma_theta^2)));
s3=s1' * s2;
p_pre_smxy = s3;
s1=(exp(-(bsxfun(@minus, edges{1}, taupost_tot)).^2/(2 * sigma_t^2)));
s2=(exp(-(bsxfun(@minus, edges{2}, thetapost_tot(Ipost)).^2)/(2*sigma_theta^2)));
s3=s1' * s2;
p_post_smxy = s3;

normtype = 'shuffstd'

if strcmp(normtype, 'unnormalised')


figure 
cmax = 0.05;
subplot(3,1,1)
imagesc(edges{1}, edges{2} * 180 / pi, p_pre_smxy')
h = colorbar;
%caxis([0,cmax])
axis([0,1250,76,90])
subplot(3,1,2)
imagesc(edges{1}, edges{2} * 180 / pi,p_post_smxy')
h = colorbar;
%caxis([0,cmax])
axis([0,1250,76,90])
subplot(3,1,3)
imagesc(edges{1}, edges{2} * 180 / pi,(p_post_smxy' - p_pre_smxy'))
h = colorbar;
%caxis([-cmax,cmax])
axis([0,1250,76,90])

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



cmax1 = 10;
cmax2 = 3;

figure 
subplot(3,1,1)
imagesc(edges{1}, edges{2} * 180 / pi, p_pre_smxy' ./ normvar_pre')
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
xlabel('Time constant (ms)')
ylabel('Angle from optimal input SNR (deg)')
h = colorbar;
title(h, 'Normalized count')
caxis([0,cmax1])
axis([0,1250,76,90])
subplot(3,1,2)
imagesc(edges{1}, edges{2} * 180 / pi,p_post_smxy' ./ normvar_post')
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
xlabel('Time constant (ms)')
ylabel('Angle from optimal input SNR (deg)')
h = colorbar;
title(h, 'Normalized count')
caxis([0,cmax1])
axis([0,1250,76,90])
subplot(3,1,3)
imagesc(edges{1}, edges{2} * 180 / pi,(p_post_smxy' - p_pre_smxy') ./ normvar_diff')
set(gca, 'ydir', 'normal')
set(gca, 'fontsize', 18)
xlabel('Time constant (ms)')
ylabel('Angle from optimal input SNR (deg)')
h = colorbar;
title(h, 'Normalized count')
caxis([-cmax2,cmax2])
axis([0,1250,76,90])
hold on
X1 = (p_post_smxy' - p_pre_smxy' > p_top_prc');
X2 = (p_post_smxy' - p_pre_smxy' < p_bottom_prc');
contour(edges{1},edges{2} * 180/pi,X1, [0,1], 'color', [0.5,0.5,0.5], 'linewidth', 4, 'linestyle', '--');
contour(edges{1},edges{2} * 180/pi,X2, [0,1], 'k', 'linewidth', 4, 'linestyle', '--');

end

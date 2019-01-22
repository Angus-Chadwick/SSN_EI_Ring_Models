
if task='learning'

elseif task='switching'
    
xpre = xrel;
xpost = xirrel;
residualpre = residualrel;
residualpost = residualirrel;

end
    
cond = 'differential'

for rix=1:length(RXS)

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
    
          SNRvec = 1;
          if SNRvec
          
            Covres_pre{rix} = cov(residualpre{RXS(rix)}');
            Covres_post{rix} = cov(residualpost{RXS(rix)}');
            Covresinv_pre{rix} = inv(Covres_pre{rix});
            Covresinv_post{rix} = inv(Covres_post{rix});

            dIn_pre{rix} = Covresinv_pre{rix} * dIn_pre{rix};
            dIn_post{rix} = Covresinv_post{rix} * dIn_post{rix};
            
          end
    
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
   
    
    I = find(imag(diag(Dpre{rix})) == 0);
    Xpre{rix} = abs(Vinv_pre{rix}(I,:) * dIn_pre{rix}) / norm(dIn_pre{rix});
    D = diag(Dpre{rix});
    Ypre{rix} = D(I);
    
    I = find(imag(diag(Dpost{rix})) == 0);
    Xpost{rix} = abs(Vinv_post{rix}(I,:) * dIn_post{rix}) / norm(dIn_post{rix});
    D = diag(Dpost{rix});
    Ypost{rix} = D(I);
    
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

    scatter(taupre_tot, 180 / pi * thetapre_tot(Ipre), 100, '+', 'markeredgecolor', 'k')
    scatter(taupost_tot, 180 / pi * thetapost_tot(Ipost), 100, '+', 'markeredgecolor', 'b')    
    
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
p_pre = hist3([taupre_tot, thetapre_tot(Ipre)], edges);
p_post = hist3([taupost_tot, thetapost_tot(Ipost)], edges);

sigma_t = 100;
sigma_theta = 1 * pi / 180;

for i=1:length(edges{1})
    for j=1:length(edges{2})
p_pre_smxy(i,j) = sum(exp(-(edges{1}(i) - taupre_tot).^2/(2 * sigma_t^2) - (edges{2}(j) - thetapre_tot(Ipre)).^2 / (2*sigma_theta^2)));
p_post_smxy(i,j) = sum(exp(-(edges{1}(i) - taupost_tot).^2/(2 * sigma_t^2) - (edges{2}(j) - thetapost_tot(Ipost)).^2 / (2*sigma_theta^2)));

    end
end


figure 
cmax = 0.05;
subplot(3,1,1)
imagesc(edges{1}, edges{2} * 180 / pi, p_pre_smxy')
h = 'colorbar';
%caxis([0,cmax])
axis([0,1250,76,90])
subplot(3,1,2)
imagesc(edges{1}, edges{2} * 180 / pi,p_post_smxy')
h = 'colorbar';
%caxis([0,cmax])
axis([0,1250,76,90])
subplot(3,1,3)
imagesc(edges{1}, edges{2} * 180 / pi,p_post_smxy' - p_pre_smxy')
h = 'colorbar';
%caxis([-cmax,cmax])
axis([0,1250,76,90])
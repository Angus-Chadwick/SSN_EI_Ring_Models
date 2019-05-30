%% Compute matrix exponentials applied to information pulses

condition = 'differential'

for rix=1:length(RXS)

    % compute signal vectors
    
          MeanVpre = mean(xpre{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9):(length(CELLLAB_ALL{RXS(rix)}) + 17)), 2);
          MeanApre = mean(xpre{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9 + 17):(length(CELLLAB_ALL{RXS(rix)}) + 17 + 17)), 2);

          MeanVpost = mean(xpost{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9):(length(CELLLAB_ALL{RXS(rix)}) + 17)), 2);
          MeanApost = mean(xpost{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9 + 17):(length(CELLLAB_ALL{RXS(rix)}) + 17 + 17)), 2);
  
    

    
          if strmatch(condition, 'differential')

            dIn_pre{rix} = MeanVpre - MeanApre;
            dIn_post{rix} = MeanVpost - MeanApost;
    
            dMpre{rix} = squeeze(nanmean(nanmean(rmat_pre{RXS(rix), 1}(:,:,rngSTM), 3))) -  squeeze(nanmean(nanmean(rmat_pre{RXS(rix), 2}(:,:,rngSTM), 3)));
            SI_Out_pre{rix} = dMpre{rix} ./ PoolVar_pre_Out{RXS(rix)};
            dMpost{rix} = squeeze(nanmean(nanmean(rmat_post{RXS(rix), 1}(:,:,rngSTM), 3))) -  squeeze(nanmean(nanmean(rmat_post{RXS(rix), 2}(:,:,rngSTM), 3)));
            SI_Out_post{rix} = dMpost{rix} ./ PoolVar_post_Out{RXS(rix)};
          
            
          elseif strmatch(condition, 'Vert') 
              
            dIn_pre{rix} = MeanVpre;
            dIn_post{rix} = MeanVpost;
            
          elseif strmatch(condition, 'Ang') 
              
            dIn_pre{rix} = MeanApre;
            dIn_post{rix} = MeanApost;
            
         
          elseif strmatch(condition, 'Sum')

            dIn_pre{rix} = (MeanVpre + MeanApre)/2;
            dIn_post{rix} = (MeanVpost + MeanApost)/2;
            
            dMpre{rix} = (squeeze(nanmean(nanmean(rmat_pre{RXS(rix), 1}(:,:,rngSTM), 3))) +  squeeze(nanmean(nanmean(rmat_pre{RXS(rix), 2}(:,:,rngSTM), 3))))/2;
            dMpost{rix} = (squeeze(nanmean(nanmean(rmat_post{RXS(rix), 1}(:,:,rngSTM), 3))) +  squeeze(nanmean(nanmean(rmat_post{RXS(rix), 2}(:,:,rngSTM), 3))))/2;
          
            
          end
    
            
            Covres_pre{rix} = cov(residualpre{RXS(rix)}');
            Covres_post{rix} = cov(residualpost{RXS(rix)}');
            Covresinv_pre{rix} = inv(Covres_pre{rix});
            Covresinv_post{rix} = inv(Covres_post{rix});
            
            InputDiscriminant_pre{rix} = Covresinv_pre{rix} * dIn_pre{rix};
            InputDiscriminant_post{rix} = Covresinv_post{rix} * dIn_post{rix};
            
            ILC_In_pre(rix) = dIn_pre{rix}' * Covres_pre{rix} * dIn_pre{rix} / (norm(dIn_pre{rix})^2);
            ILC_In_post(rix) = dIn_post{rix}' * Covres_post{rix} * dIn_post{rix} / (norm(dIn_post{rix})^2);

            
            

%% compute response covariance

clear Rpre Rpost Rpre_meansub Rpost_meansub Rpre_mean Rpsot_mean 

for stim=1:3
   Rpre{stim}  = rmat_pre{RXS(rix), stim}(:,:, rngSTM);
   Rpost{stim} = rmat_post{RXS(rix), stim}(:,:, rngSTM);


   Rpre_meansub{stim}  = permute(Rpre{stim}  - repmat(nanmean(Rpre{stim}), [size(Rpre{stim},1),1,1]), [3,1,2]);
   Rpost_meansub{stim} = permute(Rpost{stim} - repmat(nanmean(Rpost{stim}), [size(Rpost{stim},1),1,1]), [3,1,2]);

   

           Rpre_mean{stim}  = reshape(Rpre_meansub{stim},  [size(Rpre{stim},1) * length(9:17), size(Rpre{stim},2)]);
           Rpost_mean{stim} = reshape(Rpost_meansub{stim}, [size(Rpost{stim},1) * length(9:17), size(Rpost{stim},2)]);

end
          

           
           Rpre_mean_tot = vertcat(Rpre_mean{:});
           Rpost_mean_tot = vertcat(Rpost_mean{:});

                      
           [Ipre, Jpre] = find(isnan(Rpre_mean_tot));
           [Ipost, Jpost] = find(isnan(Rpost_mean_tot));

           Nfrac_pre(rix) = length(unique(Ipre)) / size(Rpre_mean_tot,1);
           Nfrac_post(rix) = length(unique(Ipost)) / size(Rpost_mean_tot,1);

           Rpre_mean_tot(Ipre,:) = [];
           Rpost_mean_tot(Ipost,:) = [];

           covmat_pre{rix} = cov(Rpre_mean_tot);
           covmat_post{rix} = cov(Rpost_mean_tot);

            OutputDiscriminant_pre{rix} = inv(covmat_pre{rix}) * dMpre{rix}';
            OutputDiscriminant_post{rix} = inv(covmat_post{rix}) * dMpost{rix}';
            
            ILC_Out_pre(rix) = dMpre{rix} * covmat_pre{rix} * dMpre{rix}' / (norm(dMpre{rix})^2);
            ILC_Out_post(rix) = dMpost{rix} * covmat_post{rix} * dMpost{rix}' / (norm(dMpost{rix})^2);
            
            
    %% compute eigendecomposition of weight matrix
          
    Apre{rix} = xpre{RXS(rix)}(:, 1:size(xpre{RXS(rix)},1));
    Apost{rix} = xpost{RXS(rix)}(:, 1:size(xpost{RXS(rix)},1));
    
    [Vpre,Dpre] = eig(Apre{rix});
    [Vpost,Dpost] = eig(Apost{rix});

    Ipre  = find(real(diag(Dpre)) + 1 > 0.95);
    Ipost = find(real(diag(Dpost)) + 1 > 0.95);

        
          P_pre{rix}(:,1)  =  InputDiscriminant_pre{rix};
          P_post{rix}(:,1)  = InputDiscriminant_post{rix};
    
%           Lsum_pre{rix}(1) = 1;
%           Lsum_post{rix}(1) = 1;
          
          Vinv_pre = inv(Vpre);
          Vinv_post = inv(Vpost);
          
          Defectiveness_pre(rix)  = cond(Vpre) ;
          Defectiveness_post(rix) = cond(Vpost) ;
                    
          Vpre(:,Ipre) = 0;
          Vpost(:,Ipost) = 0;
          Vinv_pre(Ipre,:) = 0;
          Vinv_post(Ipost,:) = 0;
       
          
          Nt = 41;
    for i=1:Nt
        
            Evals_pre = (diag(Dpre));
            Evals_pre(Ipre) = nan;
            Evals_post = (diag(Dpost));
            Evals_post(Ipost) = nan;
% 
%           indspre  = find(~ismember(1:length(Evals_pre), Ipre));
%           indspost = find(~ismember(1:length(Evals_post), Ipost));
% 
%           maxpre  = find(real(Evals_pre(indspre)) == max(real(Evals_pre(indspre))));
%           maxpost = find(real(Evals_post(indspost))  == max(real(Evals_post(indspost))));

%           
%           
%           if numel(maxpre) == 1
%           
%           Lsum_pre{rix}(:,i)  = (Vpre(:,indspre(maxpre(1)))) * (1+(Evals_pre(indspre(maxpre(1)))))^(i-1) * (Vinv_pre(indspre(maxpre(1)),:)) * InputDiscriminant_pre{rix};
%           
%           elseif numel(maxpre) == 2 
%               
%           Lsum_pre{rix}(:,i)  = ( (Vpre(:,indspre(maxpre(1)))) * (1+(Evals_pre(indspre(maxpre(1)))))^(i-1) * (Vinv_pre(indspre(maxpre(1)),:)) + (Vpre(:,indspre(maxpre(2)))) * (1+(Evals_pre(indspre(maxpre(2)))))^(i-1) * (Vinv_pre(indspre(maxpre(2)),:)) ) * InputDiscriminant_pre{rix};
% 
%           end
%           
%           if numel(maxpost) == 1
%           
%           Lsum_post{rix}(:,i)  = (Vpost(:,indspost(maxpost(1)))) * (1+(Evals_post(indspost(maxpost(1)))))^(i-1) * (Vinv_post(indspost(maxpost(1)),:)) * InputDiscriminant_post{rix};
%           
%           elseif numel(maxpost) == 2 
%               
%           Lsum_post{rix}(:,i)  = ( (Vpost(:,indspost(maxpost(1)))) * (1+(Evals_post(indspost(maxpost(1)))))^(i-1) * (Vinv_post(indspost(maxpost(1)),:)) + (Vpost(:,indspost(maxpost(2)))) * (1+(Evals_post(indspost(maxpost(2)))))^(i-1) * (Vinv_post(indspost(maxpost(2)),:)) ) * InputDiscriminant_post{rix};
% 
%           end
%           
%           P_pre{rix}(:,i)  = Vpre * (Dpre + eye(size(Dpre)))^(i-1) * Vinv_pre * InputDiscriminant_pre{rix};
%           P_post{rix}(:,i)  = Vpost * (Dpost + eye(size(Dpost)))^(i-1) * Vinv_post * InputDiscriminant_post{rix};



          P_pre{rix}(:,i)   = Vpre  * (Dpre + eye(size(Dpre))).^(i-1)   * Vinv_pre  * InputDiscriminant_pre{rix};
          P_post{rix}(:,i)  = Vpost * (Dpost + eye(size(Dpost))).^(i-1) * Vinv_post * InputDiscriminant_post{rix};

%           OpNorm_pre(rix,i) = max(svd(Vpre * (Dpre + eye(size(Dpre))).^(i-1) * Vinv_pre)); 
%           SpecRad_pre(rix,i) = max(1+real(diag(Dpre)))^(i-1); 
%           OpNorm_post(rix,i) = max(svd(Vpost * (Dpost + eye(size(Dpost))).^(i-1) * Vinv_post)); 
%           SpecRad_post(rix,i) = max(1+real(diag(Dpost)))^(i-1); 
%           
%           FirstFactor_pre{rix}(:,i) = (Apre{rix} + eye(size(Dpre)))^(i-1) * InputDiscriminant_pre{rix} * (- max(real(Evals_pre)) + 1)^(i-1);
%           FirstFactor_post{rix}(:,i) = (Apost{rix} + eye(size(Dpost)))^(i-1) * InputDiscriminant_post{rix} * (- max(real(Evals_post)) + 1)^(i-1);

          
    end
    
    P_pre{rix}(:,1) = InputDiscriminant_pre{rix};
    P_post{rix}(:,1) = InputDiscriminant_post{rix};
    
end

for rix=1:length(RXS)
subplot(4,2,rix)

Xpre = real(P_pre{rix}') * OutputDiscriminant_pre{rix};
plot([0:0.125:((Nt-1)*0.125)], real(P_pre{rix}') * OutputDiscriminant_pre{rix} ./ Xpre(1) )

hold on

Xpost = real(P_post{rix}') * OutputDiscriminant_post{rix};
plot([0:0.125:((Nt-1)*0.125)], real(P_post{rix}') * OutputDiscriminant_post{rix} ./ Xpost(1) )

axis([0,5,-0.1,1.1])
set(gca, 'fontsize', 18)

xlabel('Time (s)')
ylabel('Response')
legend('Pre', 'Post')

TempIntPre(rix) = mean(Xpre/Xpre(1));
TempIntPost(rix) = mean(Xpost/Xpost(1));

% 
% for i=1:100
% 
% [FOpre, Gpre] = fit([0:0.125:((Nt-1)*0.125)]', real(P_pre{rix}') * OutputDiscriminant_pre{rix}, 'exp1');
% [FOpost, Gpost] = fit([0:0.125:((Nt-1)*0.125)]', real(P_post{rix}') * OutputDiscriminant_post{rix}, 'exp1');
% 
% tau_pre(rix,i) = -1./FOpre.b;
% R2_pre(rix,i) = Gpre.rsquare;
% 
% tau_post(rix,i) = -1./FOpost.b;
% R2_post(rix,i) = Gpost.rsquare;
% 
% end

end

%% Plot output mode clustering

for rix=1:length(RXS)
        
[Vpost,Dpost] = eig(Apost{rix});
[Vpre,Dpre] = eig(Apre{rix});
Ipre = find(imag(diag(Dpre)) == 0);
Ipost = find(imag(diag(Dpost)) == 0);
Vinv_pre = inv(Vpre);
Vinv_post = inv(Vpost);

Template_pre  = (Vinv_pre * InputDiscriminant_pre{rix}) .* (Vpre' * OutputDiscriminant_pre{rix}) / (norm(InputDiscriminant_pre{rix}) * norm(OutputDiscriminant_pre{rix}));
Template_post = (Vinv_post * InputDiscriminant_post{rix}) .* (Vpost' * OutputDiscriminant_post{rix})  / (norm(InputDiscriminant_post{rix}) * norm(OutputDiscriminant_post{rix}));

delta_pre = diag(Dpre) - diag(Dpre)';
delta_post = diag(Dpost) - diag(Dpost)';

% delta_pre(delta_pre >= 0) = nan;
% delta_post(delta_post >= 0) = nan;

for i=1:10
a_pre(i,rix) = sum(Template_pre .* delta_pre(:,1).^(i-1) ) / Dpre(1,1)^(i-1);  % compare coefficients in taylor series (envelope vs transient growth)
a_post(i,rix) = sum(Template_post .* delta_post(:,1).^(i-1) )  / Dpost(1,1)^(i-1);
end

a_pre(:,rix) = a_pre(:,rix) / a_pre(1,rix);
a_post(:,rix) = a_post(:,rix) / a_post(1,rix);


for i=1:size(delta_pre)
    for j=1:size(delta_pre)
       FF_pre{rix}(i,j) = delta_pre(i,j) * Template_pre(j);
    end
end

for i=1:size(delta_post)
    for j=1:size(delta_post)
       FF_post{rix}(i,j) =  delta_post(i,j) * Template_post(j);
    end
end

sum(FF_pre{rix},1)

subplot(4,2,rix)
hold on
scatter(FF_pre{rix}(:), delta_pre(:))
scatter(FF_post{rix}(:), delta_post(:))

end



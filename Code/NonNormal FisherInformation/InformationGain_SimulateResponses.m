%% Compute matrix exponentials applied to information pulses


condition = 'differential';

for rix=1:length(RXS)

    %% compute signal vectors
    
    Vin_pre = xpre{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9):(length(CELLLAB_ALL{RXS(rix)}) + 17));
    Ain_pre = xpre{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9 + 17):(length(CELLLAB_ALL{RXS(rix)}) + 17 + 17));
    Vin_post = xpost{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9):(length(CELLLAB_ALL{RXS(rix)}) + 17));
    Ain_post = xpost{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9 + 17):(length(CELLLAB_ALL{RXS(rix)}) + 17 + 17));

          MeanVpre = mean(Vin_pre, 2);
          MeanApre = mean(Ain_pre, 2);

          MeanVpost = mean(Vin_post, 2);
          MeanApost = mean(Ain_post, 2);
  
    
          if strmatch(condition, 'differential')

            dIn_pre{rix} = MeanVpre - MeanApre;
            dIn_post{rix} = MeanVpost - MeanApost;
    
            dMpre{rix} = squeeze(nanmean(nanmean(rmat_pre{RXS(rix), 1}(:,:,rngSTM), 3))) -  squeeze(nanmean(nanmean(rmat_pre{RXS(rix), 2}(:,:,rngSTM), 3)));
            SI_Out_pre{rix} = dMpre{rix} ./ PoolVar_pre_Out{RXS(rix)};
            dMpost{rix} = squeeze(nanmean(nanmean(rmat_post{RXS(rix), 1}(:,:,rngSTM), 3))) -  squeeze(nanmean(nanmean(rmat_post{RXS(rix), 2}(:,:,rngSTM), 3)));
            SI_Out_post{rix} = dMpost{rix} ./ PoolVar_post_Out{RXS(rix)};
            
         
          elseif strmatch(condition, 'Sum')

            dIn_pre{rix} = (MeanVpre + MeanApre)/2;
            dIn_post{rix} = (MeanVpost + MeanApost)/2;
            
            dMpre{rix} = (squeeze(nanmean(nanmean(rmat_pre{RXS(rix), 1}(:,:,rngSTM), 3))) +  squeeze(nanmean(nanmean(rmat_pre{RXS(rix), 2}(:,:,rngSTM), 3))))/2;
            dMpost{rix} = (squeeze(nanmean(nanmean(rmat_post{RXS(rix), 1}(:,:,rngSTM), 3))) +  squeeze(nanmean(nanmean(rmat_post{RXS(rix), 2}(:,:,rngSTM), 3))))/2;
        
          end
    
          %% compute residual covariance
          
          covtype = 'AVGall'; 
                    
          if strmatch(covtype, 'AVstim')
          
            residualspre_V = ResidualV_pre_TOT{RXS(rix)}(:,9:end,:);
            residualspre_A = ResidualA_pre_TOT{RXS(rix)}(:,9:end,:);
            residualspost_V = ResidualV_post_TOT{RXS(rix)}(:,9:end,:);
            residualspost_A = ResidualA_post_TOT{RXS(rix)}(:,9:end,:);
            clear residualspre_Vtot residualspre_Atot residualspost_Vtot residualspost_Atot
            for q=1:size(residualspre_V,1)
                residualspre_Vtot(q,:) = reshape(residualspre_V(q,:,:), 1, size(residualspre_V,2) * size(residualspre_V,3));
                residualspre_Atot(q,:) = reshape(residualspre_A(q,:,:), 1, size(residualspre_A,2) * size(residualspre_A,3));
            end
            for q=1:size(residualspost_V,1)
                residualspost_Vtot(q,:) = reshape(residualspost_V(q,:,:), 1, size(residualspost_V,2) * size(residualspost_V,3));
                residualspost_Atot(q,:) = reshape(residualspost_A(q,:,:), 1, size(residualspost_A,2) * size(residualspost_A,3));
            end
              residualspre_AV = horzcat(residualspre_Vtot, residualspre_Atot);
              residualspost_AV = horzcat(residualspost_Vtot, residualspost_Atot);

              Covres_pre{rix} = cov(residualspre_AV');
              Covres_post{rix} = cov(residualspost_AV');
              
          elseif strmatch(covtype, 'AVGall')
             Covres_pre{rix} = cov(residualpre{RXS(rix)}');
             Covres_post{rix} = cov(residualpost{RXS(rix)}');

          end
            Covresinv_pre{rix} = inv(Covres_pre{rix});
            Covresinv_post{rix} = inv(Covres_post{rix});
            
            InputDiscriminant_pre{rix} = Covresinv_pre{rix} * dIn_pre{rix};
            InputDiscriminant_post{rix} = Covresinv_post{rix} * dIn_post{rix};
            


%% compute response covariance

clear Rpre Rpost Rpre_meansub Rpost_meansub Rpre_mean Rpost_mean 

covtype = 'AVstim';
% if strmatch(condition, 'V_vs_G')
%     covtype = 'VGstim';
% end

if strmatch(covtype, 'AVstim')

    for stim=1:2
       Rpre{stim}  = rmat_pre{RXS(rix), stim}(:,:, rngSTM);
       Rpost{stim} = rmat_post{RXS(rix), stim}(:,:, rngSTM);

       Rpre_meansub{stim}  = permute(Rpre{stim}  - repmat(nanmean(Rpre{stim}), [size(Rpre{stim},1),1,1]), [3,1,2]);
       Rpost_meansub{stim} = permute(Rpost{stim} - repmat(nanmean(Rpost{stim}), [size(Rpost{stim},1),1,1]), [3,1,2]);

       Rpre_mean{stim}  = reshape(Rpre_meansub{stim},  [size(Rpre{stim},1) * length(9:17), size(Rpre{stim},2)]);
       Rpost_mean{stim} = reshape(Rpost_meansub{stim}, [size(Rpost{stim},1) * length(9:17), size(Rpost{stim},2)]);
    end
    
elseif  strmatch(covtype, 'AVGall')
    
    for stim=1:3
       Rpre{stim}  = rmat_pre{RXS(rix), stim}(:,:, [rngBSL,rngSTM]);
       Rpost{stim} = rmat_post{RXS(rix), stim}(:,:, [rngBSL,rngSTM]);

       Rpre_meansub{stim}  = permute(Rpre{stim}  - repmat(nanmean(Rpre{stim}), [size(Rpre{stim},1),1,1]), [3,1,2]);
       Rpost_meansub{stim} = permute(Rpost{stim} - repmat(nanmean(Rpost{stim}), [size(Rpost{stim},1),1,1]), [3,1,2]);

       Rpre_mean{stim}  = reshape(Rpre_meansub{stim},  [size(Rpre{stim},1) * length([rngBSL,rngSTM]), size(Rpre{stim},2)]);
       Rpost_mean{stim} = reshape(Rpost_meansub{stim}, [size(Rpost{stim},1) * length([rngBSL,rngSTM]), size(Rpost{stim},2)]);
    end
elseif  strmatch(covtype, 'VGstim')
    
    for stim=[1,3]
       Rpre{stim}  = rmat_pre{RXS(rix), stim}(:,:, rngSTM);
       Rpost{stim} = rmat_post{RXS(rix), stim}(:,:, rngSTM);

       Rpre_meansub{stim}  = permute(Rpre{stim}  - repmat(nanmean(Rpre{stim}), [size(Rpre{stim},1),1,1]), [3,1,2]);
       Rpost_meansub{stim} = permute(Rpost{stim} - repmat(nanmean(Rpost{stim}), [size(Rpost{stim},1),1,1]), [3,1,2]);

       Rpre_mean{stim}  = reshape(Rpre_meansub{stim},  [size(Rpre{stim},1) * length(9:17), size(Rpre{stim},2)]);
       Rpost_mean{stim} = reshape(Rpost_meansub{stim}, [size(Rpost{stim},1) * length(9:17), size(Rpost{stim},2)]);
    end    
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
            
                        
    %% compute eigendecomposition of weight matrix
          
    Apre{rix} = xpre{RXS(rix)}(:, 1:size(xpre{RXS(rix)},1));
    Apost{rix} = xpost{RXS(rix)}(:, 1:size(xpost{RXS(rix)},1));
    
    [Vpre,Dpre] = eig(Apre{rix});
    [Vpost,Dpost] = eig(Apost{rix});

     CutOff = 0.95;
    % CutOff = Inf;
    
    Ipre  = find(real(diag(Dpre)) + 1 > CutOff);
    Ipost = find(real(diag(Dpost)) + 1 > CutOff);


          Vinv_pre = inv(Vpre);
          Vinv_post = inv(Vpost);
          

          Vpre(:,Ipre) = 0;
          Vpost(:,Ipost) = 0;
          Vinv_pre(Ipre,:) = 0;
          Vinv_post(Ipost,:) = 0;
 
          cov_init_pre = nancov(rmat_pre{RXS(rix)}(:,:,16));
          cov_init_post = nancov(rmat_post{RXS(rix)}(:,:,16));
          
          Nt = 9;

          
          prop_pre{rix}(:,1) = real(Vpre  * (Dpre + eye(size(Dpre))).^1 * Vinv_pre) * (nanmean(rmat_pre{RXS(rix),1}(:,:,rngSTM(1)-1),1) - nanmean(rmat_pre{RXS(rix),2}(:,:,rngSTM(1)-1),1))';
          prop_post{rix}(:,1) = real(Vpost  * (Dpost + eye(size(Dpost))).^1 * Vinv_post) * (nanmean(rmat_post{RXS(rix),1}(:,:,rngSTM(1)-1),1) - nanmean(rmat_post{RXS(rix),2}(:,:,rngSTM(1)-1),1))';
         
    for i=1:Nt
        
            Evals_pre = (diag(Dpre));
            Evals_pre(Ipre) = nan;
            Evals_post = (diag(Dpost));
            Evals_post(Ipost) = nan;

            Wpow_pre = Vpre  * (Dpre + eye(size(Dpre))).^(i-1)   * Vinv_pre;
            Wpow_post = Vpost * (Dpost + eye(size(Dpost))).^(i-1) * Vinv_post;
            
            P_pre{rix}(:,i)   = real(Wpow_pre)  * (Vin_pre(:,Nt-i+1) - Ain_pre(:,Nt-i+1));
            P_post{rix}(:,i)  = real(Wpow_post) * (Vin_post(:,Nt-i+1) - Ain_post(:,Nt-i+1));
            
            prop_pre{rix}(:,i+1) = real(Vpre  * (Dpre + eye(size(Dpre))) * Vinv_pre) * prop_pre{rix}(:,i);
            prop_post{rix}(:,i+1) = real(Vpost  * (Dpost + eye(size(Dpost))) * Vinv_post) * prop_post{rix}(:,i);
                        
    end
    
     P_pre{rix}(:,1) = Vin_pre(:,end) - Ain_pre(:,end);
     P_post{rix}(:,1) = Vin_post(:,end) - Ain_post(:,end);

end

% input information

for rix=1:length(RXS)
InfIn_pre(rix) = (InputDiscriminant_pre{rix}' * dIn_pre{rix});
InfIn_post(rix) = (InputDiscriminant_post{rix}' * dIn_post{rix});
end

[h,p] = ttest(InfIn_pre, InfIn_post, 'tail', 'left');


% temporal integration

for rix=1:length(RXS)
Xpre(:,rix) = real(P_pre{rix}') * OutputDiscriminant_pre{rix};
Xpost(:,rix) = real(P_post{rix}') * OutputDiscriminant_post{rix};
TempIntPre(rix) = (sum(Xpre(:,rix),1) * 0.125).^2 ./ (sum(Xpre(:,rix).^2,1) * 0.125);
TempIntPost(rix) = (sum(Xpost(:,rix),1) * 0.125).^2 ./ (sum(Xpost(:,rix).^2,1) * 0.125);
end

clear InfOut_pre InfOut_post
for rix=1:length(RXS)
    
        dResp_pre{rix} = cumsum(P_pre{rix},2) + prop_pre{rix}(:,1:9);
        dResp_post{rix} = cumsum(P_post{rix},2) + prop_post{rix}(:,1:9);

 InfOut_pre(rix) = mean(dResp_pre{rix},2)' * inv(covmat_pre{rix}) * mean(dResp_pre{rix},2);
 InfOut_post(rix) = mean(dResp_post{rix},2)' * inv(covmat_post{rix}) * mean(dResp_post{rix},2);


end

InfGain_pre = (InfOut_pre ./ InfIn_pre - 1) * 100;
InfGain_post = (InfOut_post ./ InfIn_post - 1) * 100;

[h,p] = signtest(InfGain_pre, InfGain_post, 'tail', 'left')

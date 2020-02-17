%% Compute matrix exponentials applied to information pulses

condition = 'differential';

for rix=1:length(RXS)

    %% compute signal vectors
    
          MeanVpre = mean(xpre{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9):(length(CELLLAB_ALL{RXS(rix)}) + 17)), 2);
          MeanApre = mean(xpre{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9 + 17):(length(CELLLAB_ALL{RXS(rix)}) + 17 + 17)), 2);
          MeanGpre = mean(xpre{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9 + 17*2):(length(CELLLAB_ALL{RXS(rix)}) + 17*2 + 17)), 2);

          MeanVpost = mean(xpost{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9):(length(CELLLAB_ALL{RXS(rix)}) + 17)), 2);
          MeanApost = mean(xpost{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9 + 17):(length(CELLLAB_ALL{RXS(rix)}) + 17 + 17)), 2);
          MeanGpost = mean(xpost{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9 + 17*2):(length(CELLLAB_ALL{RXS(rix)}) + 17*2 + 17)), 2);
  
    
          if strmatch(condition, 'differential')

            dIn_pre{rix} = MeanVpre - MeanApre;
            dIn_post{rix} = MeanVpost - MeanApost;
           
            dMpre_time{rix} = squeeze(nanmean(rmat_pre{RXS(rix), 1})) -  squeeze(nanmean(rmat_pre{RXS(rix), 2}));
            dMpost_time{rix} = squeeze(nanmean(rmat_post{RXS(rix), 1})) -  squeeze(nanmean(rmat_post{RXS(rix), 2}));
            
            dMpre{rix} = squeeze(nanmean(nanmean(rmat_pre{RXS(rix), 1}(:,:,rngSTM), 3))) -  squeeze(nanmean(nanmean(rmat_pre{RXS(rix), 2}(:,:,rngSTM), 3)));
            SI_Out_pre{rix} = dMpre{rix} ./ PoolVar_pre_Out{RXS(rix)};
            dMpost{rix} = squeeze(nanmean(nanmean(rmat_post{RXS(rix), 1}(:,:,rngSTM), 3))) -  squeeze(nanmean(nanmean(rmat_post{RXS(rix), 2}(:,:,rngSTM), 3)));
            SI_Out_post{rix} = dMpost{rix} ./ PoolVar_post_Out{RXS(rix)};
            
         
          elseif strmatch(condition, 'Sum')

            dIn_pre{rix} = (MeanVpre + MeanApre)/2;
            dIn_post{rix} = (MeanVpost + MeanApost)/2;
            
            dMpre{rix} = (squeeze(nanmean(nanmean(rmat_pre{RXS(rix), 1}(:,:,rngSTM), 3))) +  squeeze(nanmean(nanmean(rmat_pre{RXS(rix), 2}(:,:,rngSTM), 3))))/2;
            dMpost{rix} = (squeeze(nanmean(nanmean(rmat_post{RXS(rix), 1}(:,:,rngSTM), 3))) +  squeeze(nanmean(nanmean(rmat_post{RXS(rix), 2}(:,:,rngSTM), 3))))/2;
        
          elseif strmatch(condition, 'V_vs_G')

            dIn_pre{rix} = MeanVpre - MeanGpre;
            dIn_post{rix} = MeanVpost - MeanGpost;
            
            dMpre{rix} = squeeze(nanmean(nanmean(rmat_pre{RXS(rix), 1}(:,:,rngSTM), 3)))  - squeeze(nanmean(nanmean(rmat_pre{RXS(rix), 3}(:,:,rngSTM), 3)));
            dMpost{rix} = squeeze(nanmean(nanmean(rmat_post{RXS(rix), 1}(:,:,rngSTM), 3)))  - squeeze(nanmean(nanmean(rmat_post{RXS(rix), 3}(:,:,rngSTM), 3)));                   
            
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
            
            ILC_In_pre(rix) = dIn_pre{rix}' * Covres_pre{rix} * dIn_pre{rix} / (norm(dIn_pre{rix})^2);
            ILC_In_post(rix) = dIn_post{rix}' * Covres_post{rix} * dIn_post{rix} / (norm(dIn_post{rix})^2);

            
            

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

            OutputDiscriminant_pre_time{rix} = inv(covmat_pre{rix}) * dMpre_time{rix};      
            OutputDiscriminant_post_time{rix} = inv(covmat_post{rix}) * dMpost_time{rix};      
            
            OutputDiscriminant_pre{rix} = inv(covmat_pre{rix}) * dMpre{rix}';
            OutputDiscriminant_post{rix} = inv(covmat_post{rix}) * dMpost{rix}';            
            ILC_Out_pre(rix) = dMpre{rix} * covmat_pre{rix} * dMpre{rix}' / (norm(dMpre{rix})^2);
            ILC_Out_post(rix) = dMpost{rix} * covmat_post{rix} * dMpost{rix}' / (norm(dMpost{rix})^2);
            
            
    %% compute eigendecomposition of weight matrix
          
    Apre{rix} = xpre{RXS(rix)}(:, 1:size(xpre{RXS(rix)},1));
    Apost{rix} = xpost{RXS(rix)}(:, 1:size(xpost{RXS(rix)},1));
    
    [Vpre,Dpre] = eig(Apre{rix});
    [Vpost,Dpost] = eig(Apost{rix});
    
    CutOff = 0.95;
    CutOff = Inf;

    Ipre  = find(real(diag(Dpre)) + 1 > CutOff);
    Ipost = find(real(diag(Dpost)) + 1 > CutOff);

        
          P_pre{rix}(:,1)  =  InputDiscriminant_pre{rix};
          P_post{rix}(:,1)  = InputDiscriminant_post{rix};

          
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


          P_pre{rix}(:,i)   = Vpre  * (Dpre + eye(size(Dpre))).^(i-1)   * Vinv_pre  * InputDiscriminant_pre{rix};
          P_post{rix}(:,i)  = Vpost * (Dpost + eye(size(Dpost))).^(i-1) * Vinv_post * InputDiscriminant_post{rix};

          
    end
    
    P_pre{rix}(:,1) = InputDiscriminant_pre{rix};
    P_post{rix}(:,1) = InputDiscriminant_post{rix};
    
end

%output information

for rix=1:length(RXS)
InfOut_pre(rix) = (OutputDiscriminant_pre{rix}' * dMpre{rix}');
InfOut_post(rix) = (OutputDiscriminant_post{rix}' * dMpost{rix}');
end

[h,p] = ttest(InfOut_pre, InfOut_post, 'tail', 'left')


% input information

for rix=1:length(RXS)
InfIn_pre(rix) = (InputDiscriminant_pre{rix}' * dIn_pre{rix});
InfIn_post(rix) = (InputDiscriminant_post{rix}' * dIn_post{rix});
end

[h,p] = ttest(InfIn_pre, InfIn_post, 'tail', 'left')


% temporal integration

for rix=1:length(RXS)
Xpre(:,rix) = real(P_pre{rix}') * OutputDiscriminant_pre{rix};
Xpost(:,rix) = real(P_post{rix}') * OutputDiscriminant_post{rix};
TempIntPre(rix) = (sum(Xpre(:,rix),1) * 0.125).^2 ./ (sum(Xpre(:,rix).^2,1) * 0.125);
TempIntPost(rix) = (sum(Xpost(:,rix),1) * 0.125).^2 ./ (sum(Xpost(:,rix).^2,1) * 0.125);
end


% time-varying information

for rix=1:length(RXS)
    
    OutputDiscriminant_pre{rix} = OutputDiscriminant_pre{rix} / norm(OutputDiscriminant_pre{rix});
    OutputDiscriminant_post{rix} = OutputDiscriminant_post{rix} / norm(OutputDiscriminant_post{rix});
    
          
    for i=1:size(rmat_pre{RXS(rix),1},1)
        Proj_pre{rix,1}(:,i) = OutputDiscriminant_pre{rix}' * squeeze(rmat_pre{RXS(rix),1}(i,:,:));
    end
    for i=1:size(rmat_post{RXS(rix),1},1)
        Proj_post{rix,1}(:,i) = OutputDiscriminant_post{rix}' * squeeze(rmat_post{RXS(rix),1}(i,:,:));
    end     
    for i=1:size(rmat_pre{RXS(rix),2},1)
        Proj_pre{rix,2}(:,i) = OutputDiscriminant_pre{rix}' * squeeze(rmat_pre{RXS(rix),2}(i,:,:));
    end
    for i=1:size(rmat_post{RXS(rix),2},1)
        Proj_post{rix,2}(:,i) = OutputDiscriminant_post{rix}' * squeeze(rmat_post{RXS(rix),2}(i,:,:));
    end
    
    MeanPre{rix} = nanmean(Proj_pre{rix,1},2) - nanmean(Proj_pre{rix,2},2);
    stdPre{rix} = sqrt(0.5*(nanvar(Proj_pre{rix,1},[],2) + nanvar(Proj_pre{rix,2},[],2)));
    MeanPost{rix} = nanmean(Proj_post{rix,1},2) - nanmean(Proj_post{rix,2},2);
    stdPost{rix} = sqrt(0.5*(nanvar(Proj_post{rix,1},[],2) + nanvar(Proj_post{rix,2},[],2)));

    SNRpre{rix} = abs(MeanPre{rix}) ./ stdPre{rix};
    SNRpost{rix} = abs(MeanPost{rix}) ./ stdPost{rix};
    
    for i=1:size(rmat_pre{RXS(rix),1},1)
        Proj_pre_time{rix,1}(:,i) = diag(OutputDiscriminant_pre_time{rix}' * squeeze(rmat_pre{RXS(rix),1}(i,:,:)));
    end
    for i=1:size(rmat_post{RXS(rix),1},1)
        Proj_post_time{rix,1}(:,i) = diag(OutputDiscriminant_post_time{rix}' * squeeze(rmat_post{RXS(rix),1}(i,:,:)));
    end     
    for i=1:size(rmat_pre{RXS(rix),2},1)
        Proj_pre_time{rix,2}(:,i) = diag(OutputDiscriminant_pre_time{rix}' * squeeze(rmat_pre{RXS(rix),2}(i,:,:)));
    end
    for i=1:size(rmat_post{RXS(rix),2},1)
        Proj_post_time{rix,2}(:,i) = diag(OutputDiscriminant_post_time{rix}' * squeeze(rmat_post{RXS(rix),2}(i,:,:)));
    end
    
    
    MeanPre_time{rix} = nanmean(Proj_pre_time{rix,1},2) - nanmean(Proj_pre_time{rix,2},2);
    stdPre_time{rix} = sqrt(0.5*(nanvar(Proj_pre_time{rix,1},[],2) + nanvar(Proj_pre_time{rix,2},[],2)));
    MeanPost_time{rix} = nanmean(Proj_post_time{rix,1},2) - nanmean(Proj_post_time{rix,2},2);
    stdPost_time{rix} = sqrt(0.5*(nanvar(Proj_post_time{rix,1},[],2) + nanvar(Proj_post_time{rix,2},[],2)));

    SNRpre_time{rix} = abs(MeanPre_time{rix}) ./ stdPre_time{rix};
    SNRpost_time{rix} = abs(MeanPost_time{rix}) ./ stdPost_time{rix};


end

MeanPre_meantot = mean(horzcat(MeanPre{:}),2);
MeanPost_meantot = mean(horzcat(MeanPost{:}),2);
stdPre_meantot = mean(horzcat(stdPre{:}),2);
stdPost_meantot = mean(horzcat(stdPost{:}),2);

MeanPre_SEMtot = std(horzcat(MeanPre{:}),[],2) / sqrt(length(RXS));
MeanPost_SEMtot = std(horzcat(MeanPost{:}),[],2) / sqrt(length(RXS));
stdPre_SEMtot = std(horzcat(stdPre{:}),[],2) / sqrt(length(RXS));
stdPost_SEMtot = std(horzcat(stdPost{:}),[],2) / sqrt(length(RXS));

SNRpre_meantot = mean(horzcat(SNRpre{:}),2);
SNRpost_meantot = mean(horzcat(SNRpost{:}),2);
SNRpre_SEMtot = std(horzcat(SNRpre{:}),[],2) / sqrt(length(RXS));
SNRpost_SEMtot = std(horzcat(SNRpost{:}),[],2) / sqrt(length(RXS));


subplot(3,1,1)
hold on

shadedErrorBar(tsamples([rngBSL, rngSTM]+1), SNRpre_meantot([rngBSL, rngSTM]), SNRpre_SEMtot([rngBSL, rngSTM]), 'lineprops', 'b')
shadedErrorBar(tsamples([rngBSL, rngSTM]+1), SNRpost_meantot([rngBSL, rngSTM]), SNRpost_SEMtot([rngBSL, rngSTM]), 'lineprops', 'r')
set(gca, 'fontsize', 18)
xlabel('Time from stim onset (s)')
ylabel('Response SNR')


subplot(3,1,2)
hold on

shadedErrorBar(tsamples([rngBSL, rngSTM]+1), MeanPre_meantot([rngBSL, rngSTM]), MeanPre_SEMtot([rngBSL, rngSTM]), 'lineprops', 'b')
shadedErrorBar(tsamples([rngBSL, rngSTM]+1), MeanPost_meantot([rngBSL, rngSTM]), MeanPost_SEMtot([rngBSL, rngSTM]), 'lineprops', 'r')
set(gca, 'fontsize', 18)
xlabel('Time from stim onset (s)')
ylabel('Response signal')

subplot(3,1,3)
hold on
shadedErrorBar(tsamples([rngBSL, rngSTM]+1), stdPre_meantot([rngBSL, rngSTM]), stdPre_SEMtot([rngBSL, rngSTM]), 'lineprops', 'b')
shadedErrorBar(tsamples([rngBSL, rngSTM]+1), stdPost_meantot([rngBSL, rngSTM]), stdPost_SEMtot([rngBSL, rngSTM]), 'lineprops', 'r')
set(gca, 'fontsize', 18)
xlabel('Time from stim onset (s)')
ylabel('Response noise')

%% Check for changes in preferred and anti-preferred stims




for rix=1:length(RXS)

    for c=1:size(rmat_pre{RXS(rix)}, 2)
    
        dV_sig{rix}(c) = ranksum(squeeze(nanmean(rmat_pre{RXS(rix), 1}(:,c,rngSTM), 3)),  squeeze(nanmean(rmat_post{RXS(rix), 1}(:,c,rngSTM), 3)));
        dA_sig{rix}(c) = ranksum(squeeze(nanmean(rmat_pre{RXS(rix), 2}(:,c,rngSTM), 3)),  squeeze(nanmean(rmat_post{RXS(rix), 2}(:,c,rngSTM), 3)));

        SI_Out_pre_sig{rix}(c) = ranksum(squeeze(nanmean(rmat_pre{RXS(rix), 1}(:,c,rngSTM), 3)),  squeeze(nanmean(rmat_pre{RXS(rix), 2}(:,c,rngSTM), 3)));
        SI_Out_post_sig{rix}(c) = ranksum(squeeze(nanmean(rmat_post{RXS(rix), 1}(:,c,rngSTM), 3)),  squeeze(nanmean(rmat_post{RXS(rix), 2}(:,c,rngSTM), 3)));

    end
    
          MeanVpre{rix} = squeeze(nanmean(nanmean(rmat_pre{RXS(rix), 1}(:,:,rngSTM), 3)))' ;
          MeanApre{rix} = squeeze(nanmean(nanmean(rmat_pre{RXS(rix), 2}(:,:,rngSTM), 3)))' ;

          MeanVpost{rix} = squeeze(nanmean(nanmean(rmat_post{RXS(rix), 1}(:,:,rngSTM), 3)))' ;
          MeanApost{rix} = squeeze(nanmean(nanmean(rmat_post{RXS(rix), 2}(:,:,rngSTM), 3)))' ;          

          DeltaV{rix} = MeanVpost{rix} - MeanVpre{rix};
          DeltaA{rix} = MeanApost{rix} - MeanApre{rix};
          DeltaAVpre{rix} = MeanVpre{rix} - MeanApre{rix};
          DeltaAVpost{rix} = MeanVpost{rix} - MeanApost{rix};

end   

cells = find(CELLLAB_TOT == 1);

SIpre_sigtot = horzcat(SI_Out_pre_sig{:});
SIpost_sigtot = horzcat(SI_Out_post_sig{:});
SIpre_sigtot = SIpre_sigtot(cells);
SIpost_sigtot = SIpost_sigtot(cells);

dA_sigtot = horzcat(dA_sig{:});
dV_sigtot = horzcat(dV_sig{:});

dA_sigtot = dA_sigtot(cells);
dV_sigtot = dV_sigtot(cells);

Poolvarpre_tot = horzcat( PoolVar_pre_Out{:});
Poolvarpost_tot = horzcat( PoolVar_post_Out{:});

Poolvarpre_tot = Poolvarpre_tot(cells);
Poolvarpost_tot = Poolvarpost_tot(cells);

Poolvar_tot = sqrt(0.5*(Poolvarpre_tot.^2 + Poolvarpost_tot.^2));

DeltaAVpre_tot = vertcat(DeltaAVpre{:}) ;
DeltaAVpost_tot = vertcat(DeltaAVpost{:}) ;
DeltaV_tot = vertcat(DeltaV{:})  ;
DeltaA_tot = vertcat(DeltaA{:}) ;

DeltaAVpre_tot = DeltaAVpre_tot(cells);
DeltaAVpost_tot = DeltaAVpost_tot(cells);
DeltaV_tot = DeltaV_tot( cells);
DeltaA_tot = DeltaA_tot(cells);

prefVpre = find((DeltaAVpre_tot > 0));
prefApre = find((DeltaAVpre_tot < 0));
prefVpost = find((DeltaAVpost_tot > 0));
prefApost = find((DeltaAVpost_tot < 0));

prefVpreandpost = find(and(DeltaAVpre_tot > 0, DeltaAVpost_tot > 0));
prefApreandpost = find(and(DeltaAVpre_tot < 0, DeltaAVpost_tot < 0));


Delta_Same = [DeltaV_tot(prefVpreandpost); DeltaA_tot(prefApreandpost)];
Delta_Opp = [DeltaV_tot(prefApreandpost); DeltaA_tot(prefVpreandpost)];

Delta_VV = DeltaV_tot(prefVpreandpost);
Delta_AA = DeltaA_tot(prefApreandpost);
Delta_AV = DeltaA_tot(prefVpreandpost);
Delta_VA = DeltaV_tot(prefApreandpost);

figure
subplot(1,2,1)
bar([mean(Delta_VV), mean(Delta_VA); mean(Delta_AV), mean(Delta_AA)])
set(gca, 'fontsize', 18)
set(gca, 'xtick', [1,2])
set(gca, 'xticklabel', {'Prefer V', 'Prefer A'})
legend('Stim V', 'Stim A')
ylabel('Mean response change (post - pre)')
subplot(1,2,2)
bar([var(Delta_VV), var(Delta_VA); var(Delta_AV), var(Delta_AA)])
set(gca, 'fontsize', 18)
set(gca, 'xtick', [1,2])
set(gca, 'xticklabel', {'Prefer V', 'Prefer A'})
legend('Stim V', 'Stim A')
ylabel('Var of response change (post - pre)')


subplot(2,1,1)
hold on
bar(1:2,[mean(Delta_Same), mean(Delta_Opp)])
er = errorbar(1:2, [mean(Delta_Same), mean(Delta_Opp)], [std(Delta_Same) / sqrt(length(Delta_Same)), std(Delta_Opp) / sqrt(length(Delta_Opp))], 'linewidth', 3);
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
set(gca, 'fontsize', 18)
set(gca, 'xtick', [1,2])
set(gca, 'xticklabel', {'Preferred stim', 'Non-preferred stim'})
ylabel('Mean response change (post - pre)')

subplot(2,1,2)
hold on
bar(1:2,[var(Delta_Same), var(Delta_Opp)]) % use standard error of standard deviation?
er = errorbar(1:2, [var(Delta_Same), var(Delta_Opp)], [var(Delta_Same) / sqrt(0.5*(length(Delta_Same)-1)), std(Delta_Opp) / sqrt(0.5*(length(Delta_Opp)-1))], 'linewidth', 3);
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
set(gca, 'fontsize', 18)
set(gca, 'xtick', [1,2])
set(gca, 'xticklabel', {'Preferred stim', 'Non-preferred stim'})
ylabel('Variance of response change (post - pre)')
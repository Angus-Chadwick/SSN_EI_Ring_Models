%% Equations to do regression on linear dynamical system from individual trials

clear all;

fixvars = 'Mixed';

instantaneous = 0;
instantaneous_vonly = 1;

%ROOT     = '/nfs/nhome/live/angus/Documents/Interneuron_Data/Learning/';
ROOT = '/Users/anguschadwick/Documents/Gatsby_Computer_Files/Documents/Interneuron_Data/Learning/';RSG.savedirforAngus = [ROOT,'Saved_Data/'];

if ~exist('TSK'), TSK=[]; end
if ~isfield(TSK,'LN'),
    TSK.LN=load([RSG.savedirforAngus,'LN_TLTPDPOP.mat']);
end

INCLUDE_SD = { % rsbs_sd_pop_prepro_cutoffs
    'M70_20141022_B1'	'M70_20141028_B1' % M70_20141022_B1
    'M71_20141020_B1'	'M71_20141029_B1' % M71_20141020_B1 NB small number of trials
    'M72_20141020_B1'	'M72_20141101_B1' % M72_20141020_B1
    'M73_20141020_B1'	'M73_20141028_B1' % M73_20141020_B1
    'M75_20141021_B1'	'M75_20141029_B1' % M75_20141021_B1
    'M80_20141023_B1'	'M80_20141028_B1' % M80_20141023_B1
    'M81_20141021_B1'	'M81_20141029_B1' % M81_20141021_B1 NB small number of trials
    'M87_20141021_B1'	'M87_20141028_B1' % M87_20141021_B1
    'M89_20141021_B1'	'M89_20141029_B1' % M89_20141021_B1
    'M93_20141023_B1'	'M93_20141028_B1'} % M93_20141023_B1


for rix=1:size(INCLUDE_SD,1),
    
    if rix ~= 2 & rix ~= 7
    
    clear TOT;
    for PREPST=1:2,
        clear idinf;NID = 1;
        name = INCLUDE_SD{rix,PREPST}
        ix1=strfind(name,'_B');
        
        idinf{NID}.id    = name(1:ix1-1);
        idinf{NID}.block = str2num(name(ix1+2:end));
        
        clear EXTINF;
        fname=sprintf('%s%s_B%d_extinf',RSG.savedirforAngus,idinf{NID}.id,idinf{NID}.block);
        load(fname,'EXTINF');
        nam    = EXTINF.nam;
        if strcmp(nam,'LN'), idinf{NID}.type='SD'; elseif strcmp(nam,'SW'), idinf{NID}.type='SWITCH'; end
        
        RIXSES = EXTINF.RIXSES;
        EYE    = EXTINF.EYE;
        EYEINF = EXTINF.EYEINF;
        
        clear out ADDINF;
        fname=sprintf('%s%s_B%d_dat',RSG.savedirforAngus,idinf{NID}.id,idinf{NID}.block);
        load(fname,'out','ADDINF');
          
        clear ADDTRG;
        extract_ADDTRG=0;
        if extract_ADDTRG, % rsbs_extract_data_trigger
            dbstop if error;
            ADDTRG=rsbs_extract_data_trigger(out,ADDINF,idinf,NID,EYE,EYEINF,TSK.(nam),RIXSES);
            fname=sprintf('%sADDTRG2_%s_B%d',RSG.savedirforAngus,idinf{NID}.id,idinf{NID}.block);
            save(fname,'ADDTRG','-v7.3');
        else
            fname=sprintf('%sADDTRG2_alltrials_%s_B%d',RSG.savedirforAngus,idinf{NID}.id,idinf{NID}.block);
            load(fname,'ADDTRG');
        end
        
        % save data needed for later
        TOT{PREPST}.RIXSES = RIXSES;
        TOT{PREPST}.ADDTRG = ADDTRG;        
    end % for PREPST=1:2,
   
    
    
            % find same cells pre and post
        CELL{1} = TOT{1}.ADDTRG{1}.POOLCV; % PRE
        CELL{2} = TOT{2}.ADDTRG{1}.POOLCV; % PST
        CELLUNI = unique([CELL{1};CELL{2}],'rows');
        IX=NaN(size(CELLUNI,1),2);
        for i=1:size(CELLUNI,1),
            for PREPST=1:2,
                ix=find(CELL{PREPST}(:,1)==CELLUNI(i,1)&CELL{PREPST}(:,2)==CELLUNI(i,2));
                if ~isempty(ix),
                    IX(i,PREPST)=ix;
                end
            end
        end
        IX = IX(all(~isnan(IX),2),:); % all cells with pre and post data
        
    

%% Get cell labels

 POOLCVPRE = CELL{1};
 POOLCVPOST = CELL{2};
 
        % find cell type of each cell
        CELLLABPRE=NaN(size(POOLCVPRE,1),1); % index into TL
        RIXSES = TOT{1}.RIXSES;
        for chn=1:size(POOLCVPRE,1),
            ix=find( (TSK.(nam).TL.rix==RIXSES(1))&(TSK.(nam).TL.ses==RIXSES(2))...
                &(TSK.(nam).TL.lay==POOLCVPRE(chn,1))&(TSK.(nam).TL.B==POOLCVPRE(chn,2)) );  
            if not(isempty(ix)),
                CELLLABPRE(chn) = TSK.(nam).TL.labJ(ix);
            else
                %warning(sprintf('cannot find chn %d',chn));
            end
        end

                % find cell type of each cell
        CELLLABPOST=NaN(size(POOLCVPOST,1),1); % index into TL
                RIXSES = TOT{2}.RIXSES;
        for chn=1:size(POOLCVPOST,1),
            ix=find( (TSK.(nam).TL.rix==RIXSES(1))&(TSK.(nam).TL.ses==RIXSES(2))...
                &(TSK.(nam).TL.lay==POOLCVPOST(chn,1))&(TSK.(nam).TL.B==POOLCVPOST(chn,2)) );  
            if not(isempty(ix)),
                CELLLABPOST(chn) = TSK.(nam).TL.labJ(ix);
            else
                %warning(sprintf('cannot find chn %d',chn));
            end
        end

        CELLLAB = CELLLABPRE(IX(:,1)); % common cell labels
        
        CELLLAB_ALL{rix} = CELLLAB;
                
% Inputs:
% rmat: matrix of time samples x trials x cells for the dF/F signal
% drmat: matrix of time samples x trials x cells for the change in dF/F between consecutive time samples
% Nsamples: number of time samples
% Ntrials: number of trials

% Choose behavioural conditions to include for model fitting

%Conds = 7; % grey onset only
%Conds = [1,2]; % gratings only
Conds = [1,2,7];  % Vertical, angled and grey onset
Nconds = length(Conds);


for Cond=1:Nconds

    Condi = Conds(Cond);
    
rmat_pre{rix, Cond} = TOT{1}.ADDTRG{1}.PL{1}.DCOR{Condi}.group1(:, IX(:,1), :);
rmat_post{rix, Cond} = TOT{2}.ADDTRG{1}.PL{1}.DCOR{Condi}.group1(:, IX(:,2), :);

drmat_pre{rix, Cond} = diff(rmat_pre{rix, Cond},1,3);
drmat_post{rix, Cond} = diff(rmat_post{rix, Cond},1,3);

vmatpre{rix, Cond} = TOT{1}.ADDTRG{1}.PL{1}.DCOR{Condi}.group2;
vmatpost{rix, Cond} = TOT{2}.ADDTRG{1}.PL{1}.DCOR{Condi}.group2;

% Outputs:
% Inputsmat: matrix of time samples x cells for the average input to each cell
% CONmat: matrix of cells x cells for the connectivity between cells

% Express dynamical systems equation as a matrix product

tsamples = TOT{1}.ADDTRG{1}.PL{1}.DCOR{Condi}.t;
rngSTM = TOT{1}.ADDTRG{1}.PL{1}.DCOR{Condi}.rngSTM;
rngBSL = TOT{1}.ADDTRG{1}.PL{1}.DCOR{Condi}.rngBSL;
rngtot = [rngBSL, rngSTM];
Ntrialspre{Cond} = size(drmat_pre{rix, Cond},1);
Ntrialspost{Cond} = size(drmat_post{rix, Cond},1);
Nsamples = length(rngtot);
tMOT{Cond} = TOT{1}.ADDTRG{1}.PL{1}.DCOR{Condi}.tMOT;

clear Nvpre
clear Nvpost

for Trl=1:Ntrialspre{Cond}

    vmatpre0{rix, Cond}(Trl,:)  = interp1(tMOT{Cond}, vmatpre{rix, Cond}(Trl,:), tsamples);
    Nvpre(Trl) = sum(isnan(vmatpre0{rix, Cond}(Trl,:)));

    if Nvpre(Trl) > 0
        
        V = vmatpre0{rix, Cond}(Trl,:);
        
        vmatpre0{rix, Cond}(Trl, isnan(V)) = vmatpre0{rix, Cond}(Trl, find(isnan(V)) - 1);
        
    end
    
end

for Trl=1:Ntrialspost{Cond}
    
    vmatpost0{rix, Cond}(Trl,:) = interp1(tMOT{Cond}, vmatpost{rix, Cond}(Trl,:), tsamples);
    Nvpost(Trl) = sum(isnan(vmatpost0{rix, Cond}(Trl,:)));
    
    if Nvpost(Trl) > 0
        
        V = vmatpost0{rix, Cond}(Trl,:);
        
        vmatpost0{rix, Cond}(Trl, isnan(V)) = vmatpost0{rix, Cond}(Trl, find(isnan(V)) - 1);
        
    end
        
end

if instantaneous

    drmatTOT_pre{Cond} = drmat_pre{rix, Cond}(:,:,rngtot-1); % or put -1 in this one?
    rmatTOT_pre{Cond} = rmat_pre{rix, Cond}(:,:,rngtot);
    vmatTOT_pre{Cond} = vmatpre0{rix, Cond}(:, rngtot);

    drmatTOT_post{Cond} = drmat_post{rix, Cond}(:,:,rngtot-1);
    rmatTOT_post{Cond} = rmat_post{rix, Cond}(:,:,rngtot);
    vmatTOT_post{Cond} = vmatpost0{rix, Cond}(:, rngtot);

elseif instantaneous_vonly
    
       
    drmatTOT_pre{Cond} = drmat_pre{rix, Cond}(:,:,rngtot);
    rmatTOT_pre{Cond} = rmat_pre{rix, Cond}(:,:,rngtot);
    vmatTOT_pre{Cond} = vmatpre0{rix, Cond}(:, rngtot+1);

    drmatTOT_post{Cond} = drmat_post{rix, Cond}(:,:,rngtot);
    rmatTOT_post{Cond} = rmat_post{rix, Cond}(:,:,rngtot);
    vmatTOT_post{Cond} = vmatpost0{rix, Cond}(:, rngtot+1);

    
    
else
    
    drmatTOT_pre{Cond} = drmat_pre{rix, Cond}(:,:,rngtot);
    rmatTOT_pre{Cond} = rmat_pre{rix, Cond}(:,:,rngtot);
    vmatTOT_pre{Cond} = vmatpre0{rix, Cond}(:, rngtot);

    drmatTOT_post{Cond} = drmat_post{rix, Cond}(:,:,rngtot);
    rmatTOT_post{Cond} = rmat_post{rix, Cond}(:,:,rngtot);
    vmatTOT_post{Cond} = vmatpost0{rix, Cond}(:, rngtot);

    
end
    
clear Npre
clear Npost

for i=1:Ntrialspre{Cond}

    M = squeeze(drmatTOT_pre{Cond}(i,:,:));
    Npre(i) = sum(sum(isnan(M)));

end

for i=1:Ntrialspost{Cond}

    M = squeeze(drmatTOT_post{Cond}(i,:,:));
    Npost(i) = sum(sum(isnan(M)));

end

Trialspre{Cond} = find(Npre == 0);
Trialspost{Cond} = find(Npost == 0);

drmatTOT_pre{Cond} = drmatTOT_pre{Cond}(Trialspre{Cond},:,:);
rmatTOT_pre{Cond} = rmatTOT_pre{Cond}(Trialspre{Cond},:,:);
vmatTOT_pre{Cond} = vmatTOT_pre{Cond}(Trialspre{Cond},:);

drmatTOT_post{Cond} = drmatTOT_post{Cond}(Trialspost{Cond},:,:);
rmatTOT_post{Cond} = rmatTOT_post{Cond}(Trialspost{Cond},:,:);
vmatTOT_post{Cond} = vmatTOT_post{Cond}(Trialspost{Cond},:);

Ntrialspre{Cond} = length(Trialspre{Cond});
Ntrialspost{Cond} = length(Trialspost{Cond});

end



%% Fit full (V & A, BSL & STM) simultaneously

clear rmat0pre
clear drmat0pre
clear rmat0post
clear drmat0post


for Cond = 1:Nconds

    drmatSTM_pre = drmatTOT_pre{Cond};
    rmatSTM_pre = rmatTOT_pre{Cond};
    vmatSTM_pre = vmatTOT_pre{Cond};

    drmatSTM_post = drmatTOT_post{Cond};
    rmatSTM_post = rmatTOT_post{Cond};
    vmatSTM_post = vmatTOT_post{Cond};

    rmatSTM_post_trialmean{rix, Cond} = squeeze(mean(rmatSTM_post));
    rmatSTM_pre_trialmean{rix, Cond} = squeeze(mean(rmatSTM_pre));

    M = permute(drmatSTM_pre, [3,1,2]);
    drmatSTM_pre0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]);
    M = permute(rmatSTM_pre, [3,1,2]);
    rmatSTM_pre0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series (order [trial1, trial2, trial3,...]
    M = permute(vmatSTM_pre, [2, 1]);
    vmatSTM_pre0{Cond} = reshape(M, [size(M,1) * size(M,2),1]);

    M = permute(drmatSTM_post, [3,1,2]);
    drmatSTM_post0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series
    M = permute(rmatSTM_post, [3,1,2]);
    rmatSTM_post0 = reshape(M, [size(M,1) * size(M,2), size(M,3),1]); % reshape trials and samples into long time series
    M = permute(vmatSTM_post, [2, 1]);
    vmatSTM_post0{Cond} = reshape(M, [size(M,1) * size(M,2), 1]);

    vmatpre_all{rix,Cond}   = mean(vmatTOT_pre{Cond}); 
    vmatpost_all{rix, Cond} = mean(vmatTOT_post{Cond}); 
    
    vmatpre0{rix, Cond} = vmatTOT_pre{Cond};
    vmatpost0{rix, Cond} = vmatTOT_post{Cond};

    CLS = [1,3,4,5];


     rmat0pre{Cond}= rmatSTM_pre0.';
     drmat0pre{Cond} = drmatSTM_pre0.';
 
     rmat0post{Cond} = rmatSTM_post0.';
     drmat0post{Cond} = drmatSTM_post0.';


end



%% Set up regression matrix equations

changeW = ones(4);
changeI = [1,1,1,1];

[xpre{rix}, xpost{rix}, Measurement_pre{rix}, Measurement_post{rix}, Prediction_pre{rix}, Prediction_post{rix}] = fitLDS(fixvars, changeW, changeI, drmat0pre, drmat0post, rmat0pre, rmat0post, rngtot, Ntrialspre, Ntrialspost, vmatSTM_pre0, vmatSTM_post0, CELLLAB_ALL{rix});

residualpre{rix}  = Measurement_pre{rix} - Prediction_pre{rix};
residualpost{rix} = Measurement_post{rix} - Prediction_post{rix};

%% Calculate pooled variance and selectivity of model fits


if and(ismember(1, Conds), ismember(2,Conds))  % get selectivity measures
    
    for T=1:Ntrialspre{1}

      ResidualV_pre{rix}(:,T) = mean(squeeze(drmatTOT_pre{1}(T,:,9:end)) - xpre{rix}(:, [1:size(xpre{rix},1), (size(xpre{rix},1) + 9):(size(xpre{rix},1) + 17), end]) * [squeeze(rmatTOT_pre{1}(T,:,9:end)); eye(length(9:17)); squeeze(vmatTOT_pre{1}(T,9:end))], 2);
      ResidualV_pre_TOT{rix}(:,:,T) = squeeze(drmatTOT_pre{1}(T,:,:))    - xpre{rix}(:, [1:size(xpre{rix},1), (size(xpre{rix},1) + 1):(size(xpre{rix},1) + 17), end]) * [squeeze(rmatTOT_pre{1}(T,:,:)); eye(length(1:17)); squeeze(vmatTOT_pre{1}(T,:))];
      
    end

    for T=1:Ntrialspre{2}

        ResidualA_pre{rix}(:,T) = mean(squeeze(drmatTOT_pre{2}(T,:,9:end)) - xpre{rix}(:, [1:size(xpre{rix},1), (size(xpre{rix},1) + 9 + 17):(size(xpre{rix},1) + 17 + 17), end]) * [squeeze(rmatTOT_pre{2}(T,:,9:end)); eye(length(9:17)); squeeze(vmatTOT_pre{2}(T,9:end))], 2);
        ResidualA_pre_TOT{rix}(:,:,T) = squeeze(drmatTOT_pre{2}(T,:,:))    - xpre{rix}(:, [1:size(xpre{rix},1), (size(xpre{rix},1) + 1 + 17):(size(xpre{rix},1) + 17 + 17), end]) * [squeeze(rmatTOT_pre{2}(T,:,:)); eye(length(1:17)); squeeze(vmatTOT_pre{2}(T,1:end))];

    end

    for T=1:Ntrialspost{1}

      ResidualV_post{rix}(:,T) = mean(squeeze(drmatTOT_post{1}(T,:,9:end)) - xpost{rix}(:, [1:size(xpost{rix},1), (size(xpost{rix},1) + 9):(size(xpost{rix},1) + 17), end]) * [squeeze(rmatTOT_post{1}(T,:,9:end)); eye(length(9:17)); squeeze(vmatTOT_post{1}(T,9:end))], 2);
      ResidualV_post_TOT{rix}(:,:,T) = squeeze(drmatTOT_post{1}(T,:,:))    - xpost{rix}(:, [1:size(xpost{rix},1), (size(xpost{rix},1) + 1):(size(xpost{rix},1) + 17), end]) * [squeeze(rmatTOT_post{1}(T,:,1:end)); eye(length(1:17)); squeeze(vmatTOT_post{1}(T,1:end))];
     
    end

    for T=1:Ntrialspost{2}

      ResidualA_post{rix}(:,T) = mean(squeeze(drmatTOT_post{2}(T,:,9:end)) - xpost{rix}(:, [1:size(xpost{rix},1), (size(xpost{rix},1) + 9 + 17):(size(xpost{rix},1) + 17 + 17), end]) * [squeeze(rmatTOT_post{2}(T,:,9:end)); eye(length(9:17)); squeeze(vmatTOT_post{2}(T,9:end))], 2);
      ResidualA_post_TOT{rix}(:,:,T) = squeeze(drmatTOT_post{2}(T,:,:))    - xpost{rix}(:, [1:size(xpost{rix},1), (size(xpost{rix},1) + 1 + 17):(size(xpost{rix},1) + 17 + 17), end]) * [squeeze(rmatTOT_post{2}(T,:,1:end)); eye(length(1:17)); squeeze(vmatTOT_post{2}(T,1:end))];
      
    end

      npreA = size(ResidualA_pre{rix},2);
      npreV = size(ResidualV_pre{rix},2);
      npostA = size(ResidualA_post{rix},2);
      npostV = size(ResidualV_post{rix},2);

      PoolVar_pre_Out{rix}  = sqrt(( var(mean(rmatTOT_pre{1}(:,:,9:end),3))  * (npreV-1 ) + var(mean(rmatTOT_pre{2}(:,:,9:end),3)) * (npreA-1) ) / ( (npreV - 1) + (npreA - 1)));
      PoolVar_post_Out{rix} = sqrt(( var(mean(rmatTOT_post{1}(:,:,9:end),3)) * (npostV-1) + var(mean(rmatTOT_post{2}(:,:,9:end),3)) * (npostA-1) ) / ( (npostV - 1) + (npostA - 1)));

     
end

% find selectively responding cells

    for i=1:size(rmat_pre{rix, 1}, 2)

        psel_post{rix}(i) = ranksum(nanmean(rmat_post{rix, 1}(:,i,rngtot(9:end)), 3), nanmean(rmat_post{rix,2}(:,i,rngtot(9:end)),3));
        psel_pre{rix}(i) = ranksum(nanmean(rmat_pre{rix, 1}(:,i,rngtot(9:end)), 3), nanmean(rmat_pre{rix,2}(:,i,rngtot(9:end)),3));

    end


Trialspre_all{rix} = Trialspre;
Trialspost_all{rix} = Trialspost;

drmatpre_all{rix} = horzcat(drmat0pre{:});
drmatpost_all{rix} = horzcat(drmat0post{:});

rmatpre_all{rix} = horzcat(rmat0pre{:});
rmatpost_all{rix} = horzcat(rmat0post{:});


    end
end


RXS = [1,3,4,5,6,8,9,10];
CLS = [1,3,4,5];
CELLLAB_TOT = vertcat(CELLLAB_ALL{:});


%% Generate predicted vs measured traces for vertical stimulus 
    

for rix = 1:length(RXS)
    for T=1:length(horzcat(Trialspost_all{RXS(rix)}{1}))
        
        Predicted_post{rix}(:,:,T) = Prediction_post{RXS(rix)}(:, (17 * (T-1) + 1):(17*T)) + permute(rmat_post{RXS(rix),1}(Trialspost_all{RXS(rix)}{1}(T),:, rngtot), [2,3,1]);
        Measured_post{rix}(:,:,T) = permute(rmat_post{RXS(rix),1}(Trialspost_all{RXS(rix)}{1}(T),:, rngtot+1), [2,3,1]);
        
    end
end

for rix = 1:length(RXS)
    for T=1:length(horzcat(Trialspre_all{RXS(rix)}{1}))
        
        Predicted_pre{rix}(:,:,T) = Prediction_pre{RXS(rix)}(:, (17 * (T-1) + 1):(17*T)) + permute(rmat_pre{RXS(rix),1}(Trialspre_all{RXS(rix)}{1}(T),:, rngtot), [2,3,1]);
        Measured_pre{rix}(:,:,T) = permute(rmat_pre{RXS(rix),1}(Trialspre_all{RXS(rix)}{1}(T),:, rngtot+1), [2,3,1]);
        
    end
end

figure
hold on

rix = 4;  % session number
n = 31;  % neuron number
Trl1 = 4;  % trial number
Trl2 = 30;

plot(tsamples(rngtot + 1), mean(Measured_post{rix}(n,:,:), 3), 'color', 'b', 'linewidth', 3)
plot(tsamples(rngtot + 1), Measured_post{rix}(n,:,Trl1), 'color', 'b', 'linestyle', '--', 'linewidth', 2)
plot(tsamples(rngtot + 1), Predicted_post{rix}(n,:,Trl1), 'color', 'r', 'linestyle', '--', 'linewidth', 2)

plot(tsamples(rngtot + 1), Measured_post{rix}(n,:,Trl2), 'color', 'b', 'linestyle', '-.', 'linewidth', 2)
plot(tsamples(rngtot + 1), Predicted_post{rix}(n,:,Trl2), 'color', 'r', 'linestyle', '-.', 'linewidth', 2)

box on
set(gca, 'fontsize', 22)
xlabel('Time (s)')
ylabel('dF/F')
legend('PSTH', 'Example trial 1 measured', 'Example trial 1 predicted', 'Example trial 2 measured', 'Example trial 2 predicted')

for rix = 1:length(RXS)
    for i=1:length(CELLLAB_ALL{RXS(rix)})
        
        Rsquared_pre{rix}(i) = corr(drmatpre_all{RXS(rix)}(i,:).' - residualpre{RXS(rix)}(i,:).', drmatpre_all{RXS(rix)}(i,:).').^2;
        Rsquared_post{rix}(i) = corr(drmatpost_all{RXS(rix)}(i,:).' - residualpost{RXS(rix)}(i,:).', drmatpost_all{RXS(rix)}(i,:).').^2;
        
        Rsquared_pre0{rix}(i) = 1 - sum(residualpre{RXS(rix)}(i,:).^2) / sum((drmatpre_all{RXS(rix)}(i,:) - mean(drmatpre_all{RXS(rix)}(i,:))).^2);
        Rsquared_post0{rix}(i) = 1 - sum(residualpost{RXS(rix)}(i,:).^2) / sum((drmatpost_all{RXS(rix)}(i,:) - mean(drmatpost_all{RXS(rix)}(i,:))).^2);
                
        Rsquared_pre0_r{rix}(i) = 1  - sum(residualpre{RXS(rix)}(i,:).^2) / sum((rmatpre_all{RXS(rix)}(i,:) - mean(rmatpre_all{RXS(rix)}(i,:))).^2);
        Rsquared_post0_r{rix}(i) = 1 - sum(residualpost{RXS(rix)}(i,:).^2) / sum((rmatpost_all{RXS(rix)}(i,:) - mean(rmatpost_all{RXS(rix)}(i,:))).^2);
        
    end
    
     Rsquared_pre00(rix)  = 1 - sum(sum(residualpre{RXS(rix)}.^2)) / sum(sum((drmatpre_all{RXS(rix)} - repmat(mean(drmatpre_all{RXS(rix)},2), [1, size(drmatpre_all{RXS(rix)},2)])).^2));
     Rsquared_post00(rix) = 1 - sum(sum(residualpost{RXS(rix)}.^2)) / sum(sum((drmatpost_all{RXS(rix)} - repmat(mean(drmatpost_all{RXS(rix)},2), [1, size(drmatpost_all{RXS(rix)},2)])).^2));

     Rsquared_pre00_r(rix)  = 1 - sum(sum(residualpre{RXS(rix)}.^2)) / sum(sum((rmatpre_all{RXS(rix)} - repmat(mean(rmatpre_all{RXS(rix)},2), [1, size(drmatpre_all{RXS(rix)},2)])).^2));
     Rsquared_post00_r(rix) = 1 - sum(sum(residualpost{RXS(rix)}.^2)) / sum(sum((rmatpost_all{RXS(rix)} - repmat(mean(rmatpost_all{RXS(rix)},2), [1, size(drmatpost_all{RXS(rix)},2)])).^2));

    
end

Rsquared_pre_tot = horzcat(Rsquared_pre{:});
Rsquared_post_tot = horzcat(Rsquared_post{:});


%% Relationship between weights and selectivity, with possible explanation of earlier cell-type specific prediction analyses

% Plot inputs (weighted average over sessions)

clear Inputspre_norm
clear Inputspost_norm
clear Inputspre_v_norm
clear Inputspost_v_norm
clear Meanspre
clear Meanspost
clear Meansvpre
clear Meansvpost
clear STDpre
clear STDpost
clear N


selective = 0;

CLSNAM = {'PYR', 'PV', 'SOM', 'VIP'};

for i=1:length(CLS)
    
    Inputspre{i} = [];
    Inputspost{i} = [];
    Inputspre_vel_coeff{i} = [];
    Inputspost_vel_coeff{i} = [];
    InputspreV_vel{i} = [];
    InputspostV_vel{i} = [];
    InputspreA_vel{i} = [];
    InputspostA_vel{i} = [];
             
     for rix=1:length(RXS)
         
         Inputspre{i} = vertcat(Inputspre{i}, xpre{RXS(rix)}(CELLLAB_ALL{RXS(rix)} == CLS(i), (length(CELLLAB_ALL{RXS(rix)}) + 1):(end-1)));
         Inputspost{i} = vertcat(Inputspost{i}, xpost{RXS(rix)}(CELLLAB_ALL{RXS(rix)} == CLS(i), (length(CELLLAB_ALL{RXS(rix)}) + 1):(end-1)));

         Inputspre_vel_coeff{i}  = vertcat(Inputspre_vel_coeff{i},  xpre{RXS(rix)}(CELLLAB_ALL{RXS(rix)} == CLS(i), end));
         Inputspost_vel_coeff{i} = vertcat(Inputspost_vel_coeff{i}, xpost{RXS(rix)}(CELLLAB_ALL{RXS(rix)} == CLS(i), end));
         
         VelmeanV_pre(rix) = mean(vmatpre_all{RXS(rix), 1}(9:end));
         VelmeanV_post(rix) = mean(vmatpost_all{RXS(rix), 1}(9:end));
        
         VelmeanA_pre(rix) = mean(vmatpre_all{RXS(rix), 2}(9:end));
         VelmeanA_post(rix) = mean(vmatpost_all{RXS(rix), 2}(9:end));

         InputspreV_vel{i}  = vertcat(InputspreV_vel{i},  xpre{RXS(rix)}(CELLLAB_ALL{RXS(rix)} == CLS(i), end)  * VelmeanV_pre(rix));
         InputspostV_vel{i} = vertcat(InputspostV_vel{i}, xpost{RXS(rix)}(CELLLAB_ALL{RXS(rix)} == CLS(i), end) * VelmeanV_post(rix));
         
         InputspreA_vel{i}  = vertcat(InputspreA_vel{i},  xpre{RXS(rix)}(CELLLAB_ALL{RXS(rix)} == CLS(i), end)  * VelmeanA_pre(rix));
         InputspostA_vel{i} = vertcat(InputspostA_vel{i}, xpost{RXS(rix)}(CELLLAB_ALL{RXS(rix)} == CLS(i), end) * VelmeanA_post(rix));
         
     end
                  
     if selective 
         
        Meanspre(i,:)     = mean(Inputspre{i}(psel_pre_tot(CELLLAB_TOT == CLS(i)) < 0.05, :));
        Meanspost(i,:)    = mean(Inputspost{i}(psel_post_tot(CELLLAB_TOT == CLS(i)) < 0.05, :));
        
        STDpre = std(Inputspre{i}(psel_pre_tot(CELLLAB_TOT == CLS(i)) < 0.05, :));
        STDpost = std(Inputspost{i}(psel_post_tot(CELLLAB_TOT == CLS(i)) < 0.05, :));
        SEMpre = STDpre/sqrt(sum(psel_pre_tot(CELLLAB_TOT == CLS(i)) < 0.05));
        SEMpost = STDpost/sqrt(sum(psel_post_tot(CELLLAB_TOT == CLS(i)) < 0.05));
        
     else
     
        Meanspre(i,:)     = mean(Inputspre{i});
        Meanspost(i,:)    = mean(Inputspost{i});
        
        STDpre = std(Inputspre{i});
        STDpost = std(Inputspost{i});
        SEMpre = STDpre/sqrt(size(Inputspre{i},1));
        SEMpost = STDpost/sqrt(size(Inputspost{i},1));
    
     end
       
    Meanspre_vel(i)     = mean(Inputspre_vel_coeff{i});
    Meanspost_vel(i)    = mean(Inputspost_vel_coeff{i});
    
    Nvel(i) = length(Inputspre_vel_coeff{i});
    
    STDpre_vel(i) = std(Inputspre_vel_coeff{i});
    STDpost_vel(i) = std(Inputspost_vel_coeff{i});
  
    subplot(2,2,i)
    hold on
     % errorbar(tsamples(rngtot+1), (Meanspre(i,1:17)),  SEMpre(1:17), 'color', 'b', 'linewidth', 4, 'linestyle', '--')
      errorbar(tsamples(rngtot+1), (Meanspost(i,1:17)), SEMpost(1:17),'color', 'b', 'linewidth', 4, 'linestyle', '-')
     % errorbar(tsamples(rngtot+1), (Meanspre(i,18:34)), SEMpre(18:34), 'color', 'r', 'linewidth', 4, 'linestyle', '--')
    %  errorbar(tsamples(rngtot+1), (Meanspost(i,18:34)), SEMpost(18:34),'color', 'r', 'linewidth', 4, 'linestyle', '-')
    plot(tsamples(rngtot+1), (Meanspost(i,1:17)), SEMpost(1:17),'color', 'b', 'linewidth', 4, 'linestyle', '-')
      
      plot( (tsamples(rngtot(8)) + tsamples(rngtot(9))) / 2 * [1,1], [0, 0.03], 'color', 'k', 'linewidth', 3, 'linestyle', '--')
            
   box on
   set(gca, 'fontsize', 22)
   xlabel('Time (s)')
   ylabel('External Input')
%   legend('Vertical pre', 'Vertical post', 'Angled pre', 'Angled post')
   title(CLSNAM(i))
   axis([-0.5, 1.1, 0, 0.03])
   
end
   
   figure
   hold on
   errorbar([1:4], Meanspre_vel, STDpre_vel ./ sqrt(Nvel), '.k', 'linewidth', 3)
   errorbar([1:4] + 0.1, Meanspost_vel, STDpost_vel ./ sqrt(Nvel), '.r', 'linewidth', 3)
   set(gca, 'xtick', 1:4)
   set(gca, 'xticklabel', {'PYR', 'PV', 'SOM', 'VIP'})
   box on
   set(gca, 'fontsize', 18)
   ylabel('Velocity Input Gain')
   legend('Pre', 'Post')
   
   
%% Animal by animal output selectivity

for rix=1:length(RXS)

    dMpre{rix} = squeeze(nanmean(nanmean(rmat_pre{RXS(rix), 1}(:,:,rngSTM), 3))) -  squeeze(nanmean(nanmean(rmat_pre{RXS(rix), 2}(:,:,rngSTM), 3)));
    SI_Out_pre{rix} = dMpre{rix} ./ PoolVar_pre_Out{RXS(rix)};
    
    for c=1:size(rmat_pre{RXS(rix)}, 2)
    
        SI_Out_pre_sig{rix}(c) = ranksum(squeeze(nanmean(rmat_pre{RXS(rix), 1}(:,c,rngSTM), 3)),  squeeze(nanmean(rmat_pre{RXS(rix), 2}(:,c,rngSTM), 3)));
        SI_Out_post_sig{rix}(c) = ranksum(squeeze(nanmean(rmat_post{RXS(rix), 1}(:,c,rngSTM), 3)),  squeeze(nanmean(rmat_post{RXS(rix), 2}(:,c,rngSTM), 3)));

    end
    
    dMpost{rix} = squeeze(nanmean(nanmean(rmat_post{RXS(rix), 1}(:,:,rngSTM), 3))) -  squeeze(nanmean(nanmean(rmat_post{RXS(rix), 2}(:,:,rngSTM), 3)));
    SI_Out_post{rix} = dMpost{rix} ./ PoolVar_post_Out{RXS(rix)};
    
end


%% Animal by animal input selectivity 

for rix=1:length(RXS)

          MeanVpre = mean(xpre{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9):(length(CELLLAB_ALL{RXS(rix)}) + 17)), 2);
          MeanApre = mean(xpre{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9 + 17):(length(CELLLAB_ALL{RXS(rix)}) + 17 + 17)), 2);

          MeanVpost = mean(xpost{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9):(length(CELLLAB_ALL{RXS(rix)}) + 17)), 2);
          MeanApost = mean(xpost{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9 + 17):(length(CELLLAB_ALL{RXS(rix)}) + 17 + 17)), 2);
          
          SI_In_pre{rix} = (MeanVpre - MeanApre) ./ PoolVar_pre_Out{RXS(rix)}';
          SI_In_post{rix} = (MeanVpost - MeanApost) ./ PoolVar_post_Out{RXS(rix)}';

end   

%% Calculate and Plot Mean Connectivity

CONpre = cell(4);
CONpost = CONpre;

for rix=1:length(RXS)

    Apre = xpre{RXS(rix)}(:, 1:size(xpre{RXS(rix)},1));
    Apost = xpost{RXS(rix)}(:, 1:size(xpost{RXS(rix)},1));

    Apre(1:size(Apre,1)+1:end) = NaN; % remove self-weights
    Apost(1:size(Apre,1)+1:end) = NaN; % remove self-weights
    
    for i=1:4
        for j=1:4
    
            weights = Apre(CELLLAB_ALL{RXS(rix)} == CLS(i), CELLLAB_ALL{RXS(rix)} == CLS(j));
            CONpre{i,j} = vertcat(CONpre{i,j}, weights(:));
            weights = Apost(CELLLAB_ALL{RXS(rix)} == CLS(i), CELLLAB_ALL{RXS(rix)} == CLS(j));
            CONpost{i,j} = vertcat(CONpost{i,j}, weights(:));
            
           
        end
    end
       
end

     

CONpre_Mean = cellfun(@nanmean, CONpre);
CONpost_Mean = cellfun(@nanmean, CONpost);
CONpre_SEM = cellfun(@nanstd, CONpre)./sqrt(cellfun(@length, CONpre));
CONpost_SEM = cellfun(@nanstd, CONpost)./sqrt(cellfun(@length, CONpost));


hold on
errorbar(CONpre_Mean(:), CONpre_SEM(:),'.b', 'linewidth', 4)
errorbar([1:16] + 0.1, CONpost_Mean(:), CONpost_SEM(:),'.r', 'linewidth', 4)
box on
set(gca, 'fontsize', 22)
set(gca, 'xtick', [1:16])
set(gca, 'xticklabel', {'PYR->PYR', 'PYR->PV', 'PYR->SOM', 'PYR->VIP', 'PV->PYR', 'PV->PV', 'PV->SOM', 'PV->VIP', 'SOM->PYR', 'SOM->PV', 'SOM->SOM', 'SOM->VIP', 'VIP->PYR', 'VIP->PV', 'VIP->SOM', 'VIP->VIP'})
rotateXLabels( gca, 45)
ylabel('Weight')
legend('Pre', 'Post')
title('Interaction weights before and after learning')
hold on
plot([0, 16.5], [0 0], 'color', 'k', 'linewidth', 3)
axis([0.5, 16.5, -0.005, 0.025])

    
%% Calculate noise correlations at vertical stimulus

% Test for significance of responses

NC_method = 'samples';
NC_method2 = 'jasper';

for rix=1:length(RXS)

  BSL_pre =  nanmean(rmat_pre{RXS(rix), 1}(:,:, 13:16), 3); % -0.5 to 0 s baseline
  RSP_pre =  nanmean(rmat_pre{RXS(rix), 1}(:,:, rngSTM), 3); % 0 to 1 s response

  for i=1:size(BSL_pre,2)
      
        p_resp_pre{rix}(i) = signtest(RSP_pre(:,i) - BSL_pre(:,i));
        
  end
  
  BSL_post =  nanmean(rmat_post{RXS(rix), 1}(:,:, 13:16), 3); % -0.5 to 0 s baseline
  RSP_post =  nanmean(rmat_post{RXS(rix), 1}(:,:, rngSTM), 3); % 0 to 1 s response

  for i=1:size(BSL_post,2)
      
        p_resp_post{rix}(i) = signtest(RSP_post(:,i) - BSL_post(:,i));
        
  end
  
  
end
  
corrvec_pre = cell(4);
corrvec_post = cell(4);
corrvec_pre_sig = cell(4);
corrvec_post_sig = cell(4);

for rix=1:length(RXS)
      
   Rpre  = rmat_pre{RXS(rix), 1}(:,:, rngSTM);
   Rpost = rmat_post{RXS(rix), 1}(:,:, rngSTM);

   if strmatch(NC_method,  'trialmean')
   
       Rpre = mean(Rpre, 3);
       Rpost = mean(Rpost,3);
       Rpre_mean = Rpre - repmat(nanmean(Rpre), [size(Rpre,1), 1]);
       Rpost_mean = Rpost - repmat(nanmean(Rpost), [size(Rpost,1), 1]);
   
   elseif strmatch(NC_method, 'samples')
       
       Rpre_meansub  = permute(Rpre  - repmat(nanmean(Rpre), [size(Rpre,1),1,1]), [3,1,2]);
       Rpost_meansub = permute(Rpost - repmat(nanmean(Rpost), [size(Rpost,1),1,1]), [3,1,2]);
   
       if strmatch(NC_method2, 'angus')  % my more efficient code

           Rpre_mean  = reshape(Rpre_meansub,  [size(Rpre,1) * length(tsims), size(Rpre,2)]);
           Rpost_mean = reshape(Rpost_meansub, [size(Rpost,1) * length(tsims), size(Rpost,2)]);

           [Ipre, Jpre] = find(isnan(Rpre_mean));
           [Ipost, Jpost] = find(isnan(Rpost_mean));

           Nfrac_pre(rix) = length(unique(Ipre)) / size(Rpre_mean,1);
           Nfrac_post(rix) = length(unique(Ipost)) / size(Rpost_mean,1);

           Rpre_mean(Ipre,:) = [];
           Rpost_mean(Ipost,:) = [];

           [corrmat_pre{rix}, corrmat_pre_sig{rix}] = corr(Rpre_mean);
           [corrmat_post{rix}, corrmat_post_sig{rix}] = corr(Rpost_mean);

       elseif strmatch(NC_method2, 'jasper') % a more direct replication of Jasper's code
           
           for i=1:length(CELLLAB_ALL{RXS(rix)})
               A1_pre = squeeze(Rpre_meansub(:,:,i));
               A1_post = squeeze(Rpost_meansub(:,:,i));
               for j=1:length(CELLLAB_ALL{RXS(rix)})
                   A2_pre = squeeze(Rpre_meansub(:,:,j));
                   dum = corrcov(nancov(A1_pre(:), A2_pre(:)));
                   corrmat_pre{rix}(i,j)=dum(1,2);
                   A2_post = squeeze(Rpost_meansub(:,:,j));
                   dum = corrcov(nancov(A1_post(:), A2_post(:)));
                   corrmat_post{rix}(i,j)=dum(1,2);
               end
           end
       end
           
   end
       
  
      
   corrmat_pre{rix}(p_resp_pre{rix} > 0.05, :) = NaN;
   corrmat_pre{rix}(:,p_resp_pre{rix} > 0.05) = NaN;

   corrmat_post{rix}(p_resp_post{rix} > 0.05,:) = NaN;
   corrmat_post{rix}(:, p_resp_post{rix} > 0.05) = NaN;

   
   for i=1:length(CLS)
       for j=1:length(CLS)
           
           dumpre = corrmat_pre{rix}(CELLLAB_ALL{RXS(rix)} == CLS(i),  CELLLAB_ALL{RXS(rix)} == CLS(j));
           dumpre(dumpre == 1) = nan;
           dumpost = corrmat_post{rix}(CELLLAB_ALL{RXS(rix)} == CLS(i), CELLLAB_ALL{RXS(rix)} == CLS(j));
           dumpost(dumpost == 1) = nan;
           corrvec_pre{i,j}  = vertcat(corrvec_pre{i,j},  dumpre(:));
           corrvec_post{i,j} = vertcat(corrvec_post{i,j}, dumpost(:));

           dumpre = corrmat_pre_sig{rix}(CELLLAB_ALL{RXS(rix)} == CLS(i),  CELLLAB_ALL{RXS(rix)} == CLS(j));
           dumpost = corrmat_post_sig{rix}(CELLLAB_ALL{RXS(rix)} == CLS(i), CELLLAB_ALL{RXS(rix)} == CLS(j));
           corrvec_pre_sig{i,j}  = vertcat(corrvec_pre_sig{i,j},  dumpre(:));
           corrvec_post_sig{i,j} = vertcat(corrvec_post_sig{i,j}, dumpost(:));
           
       end
   end
           
end


%% Analytical measure of total noise correlations with weight manipulations

manipulation = 'removeW1W2';

W1 = 4;
W2 = 3;
tsims = 9:17;

[ctot_pre_allweights{W1,W2}, ctot_post_allweights{W1,W2}, ctot_pre_manipW1W2{W1,W2}, ctot_post_manipW1W2{W1,W2}] = noisecorr_analytic(manipulation, W1, W2, xpre, xpost, rmat_pre, rmat_post, ResidualV_pre_TOT, ResidualV_post_TOT, vmatpre0, vmatpost0, Trialspre_all, Trialspost_all, tsims, rngtot,CELLLAB_ALL);



for i=1:length(ctot_pre_allweights)
    
    dcorr(i,1) = nanmean(ctot_pre_allweights{i} - ctot_pre_manipW1W2{i});
    dcorr(i,2) = nanmean(ctot_post_allweights{i} - ctot_post_manipW1W2{i});

    dcorr_sem(i,1) = nanstd(ctot_pre_allweights{i} - ctot_pre_manipW1W2{i}) / sqrt(length(ctot_pre_allweights{i}));
    dcorr_sem(i,2) = nanstd(ctot_post_allweights{i} - ctot_post_manipW1W2{i}) / sqrt(length(ctot_post_allweights{i}));

end
    
hold on

hb = bar(1:size(dcorr,1), -dcorr);

for ib = 1:numel(hb)

      % Find the centers of the bars
      xData = get(get(hb(ib),'Children'),'XData');
      barCenters = mean(unique(xData,'rows'));

      errorbar(barCenters,-dcorr(:,ib),dcorr_sem(:,ib),'k.', 'linewidth', 3);

end

set(gca, 'fontsize', 18)
set(gca, 'xtick', [1:10])
set(gca, 'xticklabel', {'PYR-PYR', 'VIP-VIP', 'SOM-SOM', 'PV-PV', 'PV-VIP', 'PYR-PV', 'PV-SOM', 'PYR-VIP', 'PYR-SOM', 'SOM-VIP'})
rotateXLabels( gca, 45)
box on
legend('Pre', 'Post')
axis([0.5, 10.5, -0.08, 0.04])

%% Analytical decomposition of noise correlations into sources of underlying variability

removeweights = 0;
shuffleresiduals = 0;

[cvel_pre0, cvel_post0, cres_pre0, cres_post0, cinit_pre0, cinit_post0, cvel_res_pre0, cvel_res_post0, cvel_init_pre0, cvel_init_post0, cres_init_pre0, cres_init_post0, ctot_pre, ctot_post] = noisecorr_analytic_decomp(shuffleresiduals, removeweights, xpre, xpost, rmat_pre, rmat_post, ResidualV_pre_TOT, ResidualV_post_TOT, vmatpre0, vmatpost0, Trialspre_all, Trialspost_all, tsims, rngtot, CELLLAB_ALL);

Data_Means = [cres_pre0(:), cres_post0(:), cinit_pre0(:), cinit_post0(:), cvel_pre0(:), cvel_post0(:), cvel_res_pre0(:), cvel_res_post0(:), cvel_init_pre0(:), cvel_init_post0(:), cres_init_pre0(:), cres_init_post0(:)].';
hb = bar(1:size(Data_Means,2),Data_Means');

set(gca, 'xticklabel', {'PYR-PYR', 'VIP-VIP', 'SOM-SOM', 'PV-PV', 'PV-VIP', 'PYR-PV', 'PV-SOM', 'PYR-VIP', 'PYR-SOM', 'SOM-VIP'})
set(gca, 'xtick', [1:10])
set(gca, 'fontsize', 18)
axis([0.5, 10.5, -0.04, 0.2])
legend('Residuals pre', 'Residuals post', 'Initial state pre', 'Initial state post', 'Velocity pre', 'Velocity post', 'Velocity-Residual pre', 'Velocity-Residual post', 'Velocity-Initial state pre', 'Velocity-Initial state post', 'Residual-Initial state pre', 'Residual-Initial state post')
ylabel('Contribution to total correlation')
rotateXLabels( gca, 45)

%% Analytically calculate the contribution to selectivity arising from each source

var_method = 'fullvar';
manipulation = 'removeall';

clear SI
clear poolvar

tsims = 9:17;

W1 = 1; W2 = 1;

        [SI{W1,W2}, poolvar] = selectivity_analytic(manipulation, var_method, W1, W2, xpre, xpost, rmat_pre, rmat_post, ResidualV_pre_TOT, ResidualV_post_TOT, ResidualA_pre_TOT, ResidualA_post_TOT,vmatpre0, vmatpost0, Trialspre_all, Trialspost_all, tsims, rngtot, CELLLAB_ALL);


SI_Out_pre_tot = horzcat(SI_Out_pre{:});
SI_Out_post_tot = horzcat(SI_Out_post{:});

SI_model_pre = SI{1}.SIinputs(:,1) + SI{1}.SIinitstate(:,1) + SI{1}.SIvel(:,1) ;
SI_model_post = SI{1}.SIinputs(:,2) + SI{1}.SIinitstate(:,2) + SI{1}.SIvel(:,2) ;

SI_model_pre_full = SI{1}.SIinputs(:,1) + SI{1}.SIinitstate(:,1) + SI{1}.SIvel(:,1) + SI{1}.SIres(:,1) ;
SI_model_post_full = SI{1}.SIinputs(:,2) + SI{1}.SIinitstate(:,2) + SI{1}.SIvel(:,2) + SI{1}.SIres(:,2) ;

for j=1:4
    Model_pre(j) = mean(abs(SI_model_pre(CELLLAB_TOT == CLS(j))));
    Model_post(j) = mean(abs(SI_model_post(CELLLAB_TOT == CLS(j))));
    Model_pre_full(j) = mean(abs(SI_model_pre_full(CELLLAB_TOT == CLS(j))));
    Model_post_full(j) = mean(abs(SI_model_post_full(CELLLAB_TOT == CLS(j))));
    Data_pre(j) = mean(abs(SI_Out_pre_tot(CELLLAB_TOT == CLS(j))));
    Data_post(j) = mean(abs(SI_Out_post_tot(CELLLAB_TOT == CLS(j))));
end

bar([Data_post - Data_pre; Model_post - Model_pre]')

CellNum = 3;

Xpre = [ones(size(SIinputs_rw(CELLLAB_TOT == CellNum, 1))), SIinputs_rw(CELLLAB_TOT == CellNum, 1)];
Xpost = [ones(size(SIinputs_rw(CELLLAB_TOT == CellNum, 2))), SIinputs_rw(CELLLAB_TOT == CellNum, 2)];
Ypre = SIinputs(CELLLAB_TOT == CellNum, 1);
Ypost = SIinputs(CELLLAB_TOT == CellNum, 2);
Wpre = Ypre.' / Xpre.';
Wpost = Ypost.' / Xpost.';

figure
subplot(121)
hold on
scatter(SIinputs_rw(CELLLAB_TOT == CellNum,1), SIinputs(CELLLAB_TOT == CellNum,1))
plot([-6, 6], [ -6, 6], 'k')
plot([-6,6], Wpre(2) * [-6,6] + Wpre(1))
axis([-6,6,-6,6])
set(gca, 'fontsize', 18)
title('SOM cells, pre learning')
xlabel('Selectivity of inputs without weights')
ylabel('Selectivity of inputs convolved with weights')
box on
subplot(122)
hold on
scatter(SIinputs_rw(CELLLAB_TOT == CellNum,2), SIinputs(CELLLAB_TOT == CellNum,2))
plot([-6, 6], [ -6, 6], 'k')
plot([-6,6], Wpost(2) * [-6,6] + Wpost(1))
axis([-6,6,-6,6])
set(gca, 'fontsize', 18)
title('SOM cells, post learning')
xlabel('Selectivity of inputs without weights')
ylabel('Selectivity of inputs convolved with weights')
box on

% 
%  set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20])
%  figname = strcat('Delete', CELLNAM(W2), 'to', CELLNAM(W1));
%  print('-dpng', figname{1}, '-r300'); 
% clf

%end

% Data_Means = [Mpre;Mpost;Mpre_rw;Mpost_rw];
% Data_SEM   = [Mpre_sem;Mpost_sem;Mpre_rw_sem;Mpost_rw_sem];
% 
% hb = bar(1:size(Data_Means,2),Data_Means');
% 
% % For each set of bars, find the centers of the bars, and write error bars
% 
% for ib = 1:numel(hb)
% 
%       % Find the centers of the bars
%       xData = get(get(hb(ib),'Children'),'XData');
%       barCenters = mean(unique(xData,'rows'));
% 
%       errorbar(barCenters,Data_Means(ib,:),Data_SEM(ib,:),'k.', 'linewidth', 3);
% 
% end

%% Test specificity of interaction weights between vertical and angled preferring cells



selective = 1;

for prepost = {'pre', 'post'}

    Aopp = cell(4);
    Asame = cell(4);
    
    A_AV = cell(4);
    A_VA = cell(4);
    A_AA = cell(4);
    A_VV = cell(4);
    
    for CellNum1 = 1:4
        for CellNum2 = 1:4

            for rix = 1:length(RXS)

                if strmatch(prepost, 'pre') 
                    
                    if selective
                    
                        I1 = find(CELLLAB_ALL{RXS(rix)} == CLS(CellNum1) & psel_pre{RXS(rix)}.' < 0.05);
                        I2 = find(CELLLAB_ALL{RXS(rix)} == CLS(CellNum2) & psel_pre{RXS(rix)}.' < 0.05);

                    else

                        I1 = find(CELLLAB_ALL{RXS(rix)} == CLS(CellNum1));
                        I2 = find(CELLLAB_ALL{RXS(rix)} == CLS(CellNum2));

                    end
                        
                    Apre = xpre{RXS(rix)}(:, 1:size(xpost{RXS(rix)},1));
                    Apv = Apre(I1, I2);    
                    
                    SI = SI_In_pre{rix};
                %     SI = SI_Out_pre{rix};
                
                dum = or(bsxfun(@and, SI(I1) > 0, SI(I2)' > 0), bsxfun(@and, SI(I1)<0, SI(I2)'<0));
                SamePref_pre{CellNum1, CellNum2, RXS(rix)} = zeros(length(CELLLAB_ALL{RXS(rix)}));
                SamePref_pre{CellNum1, CellNum2, RXS(rix)}(I1,I2) = dum;
                
                dum = or(bsxfun(@and, SI(I1) > 0, SI(I2)' < 0), bsxfun(@and, SI(I1)<0, SI(I2)'>0));
                OppPref_pre{CellNum1, CellNum2, RXS(rix)} = zeros(length(CELLLAB_ALL{RXS(rix)}));
                OppPref_pre{CellNum1, CellNum2, RXS(rix)}(I1,I2) = dum;
   
                elseif strmatch(prepost, 'post')
                    
                    if selective
                    
                        I1 = find(CELLLAB_ALL{RXS(rix)} == CLS(CellNum1) & psel_post{RXS(rix)}.' < 0.05);
                        I2 = find(CELLLAB_ALL{RXS(rix)} == CLS(CellNum2) & psel_post{RXS(rix)}.' < 0.05);

                    else

                        I1 = find(CELLLAB_ALL{RXS(rix)} == CLS(CellNum1));
                        I2 = find(CELLLAB_ALL{RXS(rix)} == CLS(CellNum2));

                    end
                                        
                    Apost = xpost{RXS(rix)}(:, 1:size(xpost{RXS(rix)},1));
                    Apv = Apost(I1, I2);
                    
                    SI = SI_In_post{rix};
                %    SI = SI_Out_post{rix};

                dum = or(bsxfun(@and, SI(I1) > 0, SI(I2)' > 0), bsxfun(@and, SI(I1)<0, SI(I2)'<0));
                SamePref_post{CellNum1, CellNum2, RXS(rix)} = zeros(length(CELLLAB_ALL{RXS(rix)}));
                SamePref_post{CellNum1, CellNum2, RXS(rix)}(I1,I2) = dum;
                
                dum = or(bsxfun(@and, SI(I1) > 0, SI(I2)' < 0), bsxfun(@and, SI(I1)<0, SI(I2)'>0));
                OppPref_post{CellNum1, CellNum2, RXS(rix)} = zeros(length(CELLLAB_ALL{RXS(rix)}));
                OppPref_post{CellNum1, CellNum2, RXS(rix)}(I1,I2) = dum;
   
                    
                 end

                for i=1:length(Apv)

                    Apv(i,i) = nan;

                end

                dum = Apv(SI(I1) < 0, SI(I2) > 0);
                Aopp{CellNum1, CellNum2} = vertcat(Aopp{CellNum1, CellNum2}, dum(:));
                dum = Apv(SI(I1) > 0, SI(I2) < 0);
                Aopp{CellNum1, CellNum2} = vertcat(Aopp{CellNum1, CellNum2}, dum(:));

                dum = Apv(SI(I1) < 0, SI(I2) < 0);
                Asame{CellNum1, CellNum2} = vertcat(Asame{CellNum1, CellNum2}, dum(:));
                dum = Apv(SI(I1) > 0, SI(I2) > 0);
                Asame{CellNum1, CellNum2} = vertcat(Asame{CellNum1, CellNum2}, dum(:));

                % all 4 types
                
                dum = Apv(SI(I1) < 0, SI(I2) > 0);
                A_AV{CellNum1, CellNum2} = vertcat(A_AV{CellNum1, CellNum2}, dum(:));
                dum = Apv(SI(I1) > 0, SI(I2) < 0);
                A_VA{CellNum1, CellNum2} = vertcat(A_VA{CellNum1, CellNum2}, dum(:));

                dum = Apv(SI(I1) < 0, SI(I2) < 0);
                A_AA{CellNum1, CellNum2} = vertcat(A_AA{CellNum1, CellNum2}, dum(:));
                dum = Apv(SI(I1) > 0, SI(I2) > 0);
                A_VV{CellNum1, CellNum2} = vertcat(A_VV{CellNum1, CellNum2}, dum(:));
                
            end
            
        end
    end

    if strmatch(prepost, 'pre')
        
        Asame_pre_tot = cellfun(@nanmean, Asame);
        Asame_pre_tot_sem = cellfun(@nanstd, Asame) ./ sqrt(cellfun(@length, Asame));
        Aopp_pre_tot = cellfun(@nanmean, Aopp);
        Aopp_pre_tot_sem = cellfun(@nanstd, Aopp) ./ sqrt(cellfun(@length, Aopp));
        
        A_AA_pre_tot = cellfun(@nanmean, A_AA);
        A_VV_pre_tot = cellfun(@nanmean, A_VV);
        A_AV_pre_tot = cellfun(@nanmean, A_AV);
        A_VA_pre_tot = cellfun(@nanmean, A_VA);
        
        for i=1:4
            for j=1:4
                
                  pval_pre(i,j) = ranksum(Asame{i,j}, Aopp{i,j});
                  
            end
        end
                            
    elseif strmatch(prepost, 'post')
       
        Asame_post_tot = cellfun(@nanmean, Asame);
        Asame_post_tot_sem = cellfun(@nanstd, Asame) ./ sqrt(cellfun(@length, Asame));
        Aopp_post_tot = cellfun(@nanmean, Aopp);
        Aopp_post_tot_sem = cellfun(@nanstd, Aopp) ./ sqrt(cellfun(@length, Aopp));
        
        A_AA_post_tot = cellfun(@nanmean, A_AA);
        A_VV_post_tot = cellfun(@nanmean, A_VV);
        A_AV_post_tot = cellfun(@nanmean, A_AV);
        A_VA_post_tot = cellfun(@nanmean, A_VA);
                
        for i=1:4
            for j=1:4
                
                  pval_post(i,j) = ranksum(Asame{i,j}, Aopp{i,j});
                  
            end
        end
        
    end
    
end

Data_Means = [Asame_pre_tot(:), Aopp_pre_tot(:), Asame_post_tot(:), Aopp_post_tot(:)];
Data_SEM = [Asame_pre_tot_sem(:), Aopp_pre_tot_sem(:), Asame_post_tot_sem(:), Aopp_post_tot_sem(:)];
hb = bar(1:size(Data_Means,1),Data_Means);
hold on

for ib = 1:numel(hb)
    
    xData = get(get(hb(ib),'Children'),'XData'); 
    barCenters = mean(unique(xData,'rows'));
    errorbar(barCenters,Data_Means(:,ib),Data_SEM(:,ib),'k.', 'linewidth', 2);
    
end

set(gca, 'fontsize', 18)
set(gca, 'xtick', 1:16)
set(gca, 'xticklabel', {'PYR->PYR', 'PYR->PV', 'PYR->SOM', 'PYR->VIP', 'PV->PYR', 'PV->PV', 'PV->SOM', 'PV->VIP', 'SOM->PYR', 'SOM->PV', 'SOM->SOM', 'SOM->VIP', 'VIP->PYR', 'VIP->PV', 'VIP->SOM', 'VIP->VIP'})
rotateXLabels( gca, 45)
legend('Same preference, pre','Opposite preference, pre', 'Same preference, post',  'Opposite preference, post')
ylabel('Interaction weight')


%% Delete specific weights and calculate selectivity

% delete same pref PYR to PV weights

CELL1 = 2;
CELL2 = 1;

% PYR to PV weights

for rix = RXS
    
    xpre_rwsame{rix} = xpre{rix};
    xpre_rwsame{rix}(logical(SamePref_pre{CELL1,CELL2,rix})) = 0;

    xpost_rwsame{rix} = xpost{rix};
    xpost_rwsame{rix}(logical(SamePref_post{CELL1,CELL2,rix})) = 0;
    
end

[SI_rwsame, poolvar_rwsame] = selectivity_analytic_prepostshuff(xpre_rwsame, xpost_rwsame, rmat_pre, rmat_post, ResidualV_pre_TOT, ResidualV_post_TOT, ResidualA_pre_TOT, ResidualA_post_TOT,vmatpre0, vmatpost0, Trialspre_all, Trialspost_all, tsims, rngtot, CELLLAB_ALL);


% delete opposite pref PYR to PV weights

% PYR to PV weights

for rix = RXS
    
    xpre_rwopp{rix} = xpre{rix};
    xpre_rwopp{rix}(logical(OppPref_pre{CELL1,CELL2,rix})) = 0;

    xpost_rwopp{rix} = xpost{rix};
    xpost_rwopp{rix}(logical(OppPref_post{CELL1,CELL2,rix})) = 0;
    
end

[SI_rwopp, poolvar_rwopp] = selectivity_analytic_prepostshuff(xpre_rwopp, xpost_rwopp, rmat_pre, rmat_post, ResidualV_pre_TOT, ResidualV_post_TOT, ResidualA_pre_TOT, ResidualA_post_TOT,vmatpre0, vmatpost0, Trialspre_all, Trialspost_all, tsims, rngtot, CELLLAB_ALL);


SI_rwsame_pre_tot = SI_rwsame.SIinputs(:,1) + SI_rwsame.SIinitstate(:,1) + SI_rwsame.SIvel(:,1) + SI_rwsame.SIres(:,1);
SI_rwopp_pre_tot = SI_rwopp.SIinputs(:,1) + SI_rwopp.SIinitstate(:,1) + SI_rwopp.SIvel(:,1) + SI_rwopp.SIres(:,1);
SI_rwsame_post_tot = SI_rwsame.SIinputs(:,2) + SI_rwsame.SIinitstate(:,2) + SI_rwsame.SIvel(:,2) + SI_rwsame.SIres(:,2);
SI_rwopp_post_tot = SI_rwopp.SIinputs(:,2) + SI_rwopp.SIinitstate(:,2) + SI_rwopp.SIvel(:,2) + SI_rwopp.SIres(:,2);
SI_Out_pre_tot = horzcat(SI_Out_pre{:});
SI_Out_post_tot = horzcat(SI_Out_post{:});

CellType = 3;

MEAN = [mean(abs(SI_rwopp_pre_tot(CELLLAB_TOT == CellType))) - mean(abs(SI_Out_pre_tot(CELLLAB_TOT == CellType))), mean(abs(SI_rwopp_post_tot(CELLLAB_TOT == CellType))) - mean(abs(SI_Out_post_tot(CELLLAB_TOT == CellType))); mean(abs(SI_rwsame_pre_tot(CELLLAB_TOT == CellType))) - mean(abs(SI_Out_pre_tot(CELLLAB_TOT == CellType))), mean(abs(SI_rwsame_post_tot(CELLLAB_TOT == CellType))) - mean(abs(SI_Out_post_tot(CELLLAB_TOT == CellType)))];
SEM = [std(abs(SI_rwopp_pre_tot(CELLLAB_TOT == CellType))' - abs(SI_Out_pre_tot(CELLLAB_TOT == CellType))), std(abs(SI_rwopp_post_tot(CELLLAB_TOT == CellType))' - abs(SI_Out_post_tot(CELLLAB_TOT == CellType))); std(abs(SI_rwsame_pre_tot(CELLLAB_TOT == CellType))' - abs(SI_Out_pre_tot(CELLLAB_TOT == CellType))), std(abs(SI_rwsame_post_tot(CELLLAB_TOT == CellType))' - abs(SI_Out_post_tot(CELLLAB_TOT == CellType)))] ./ sqrt(length(CELLLAB_TOT == CellType));

hold on
hb = bar(1:size(MEAN,1), MEAN);

for ib = 1:numel(hb)
    
    xData = get(get(hb(ib),'Children'),'XData'); 
    barCenters = mean(unique(xData,'rows'));
    errorbar(barCenters,MEAN(:,ib),SEM(:,ib),'k.', 'linewidth', 2);
    
end

set(gca, 'fontsize', 18)
set(gca,'xtick', 1:2)
set(gca, 'xticklabel', {'Delete Opposite', 'Delete Same'})
ylabel('\Delta abs SI')
legend('Pre Learning', 'Post Learning')
title('Effect of PYR to PV weight deletion on PV selectivity')
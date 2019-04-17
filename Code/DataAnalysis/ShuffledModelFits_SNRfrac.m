%% Equations to do regression on linear dynamical system from individual trials

clear all;

for Nshuff = 1:1000

    Nshuff

fixvars = 'Mixed';

instantaneous = 0;
instantaneous_vonly = 1;

%ROOT     = '/nfs/nhome/live/angus/Documents/Interneuron_Data/Learning/';
ROOT = '/Users/anguschadwick/Documents/Gatsby_Computer_Files/Documents/Interneuron_Data/Learning/';
RSG.savedirforAngus = [ROOT,'Saved_Data/'];

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


RandVec = randperm(Ntrialspre{Cond} + Ntrialspost{Cond});

drmatTOT{Cond} = [drmatTOT_pre{Cond}; drmatTOT_post{Cond}];
rmatTOT{Cond} = [rmatTOT_pre{Cond}; rmatTOT_post{Cond}];
vmatTOT{Cond} = [vmatTOT_pre{Cond}; vmatTOT_post{Cond}];

drmatTOT{Cond} = drmatTOT{Cond}(RandVec,:,:);
rmatTOT{Cond} = rmatTOT{Cond}(RandVec,:,:);
vmatTOT{Cond} = vmatTOT{Cond}(RandVec,:);

drmatTOT_pre{Cond} = drmatTOT{Cond}(1:Ntrialspre{Cond},:,:);
drmatTOT_post{Cond} = drmatTOT{Cond}(1:Ntrialspost{Cond},:,:);

rmatTOT_pre{Cond} = rmatTOT{Cond}(1:Ntrialspre{Cond},:,:);
rmatTOT_post{Cond} = rmatTOT{Cond}(1:Ntrialspost{Cond},:,:);

vmatTOT_pre{Cond} = vmatTOT{Cond}(1:Ntrialspre{Cond},:);
vmatTOT_post{Cond} = vmatTOT{Cond}(1:Ntrialspost{Cond},:);

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


%%

for rix=1:length(RXS)

          MeanVpre = mean(xpre{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9):(length(CELLLAB_ALL{RXS(rix)}) + 17)), 2);
          MeanApre = mean(xpre{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9 + 17):(length(CELLLAB_ALL{RXS(rix)}) + 17 + 17)), 2);

          MeanVpost = mean(xpost{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9):(length(CELLLAB_ALL{RXS(rix)}) + 17)), 2);
          MeanApost = mean(xpost{RXS(rix)}(1:length(CELLLAB_ALL{RXS(rix)}), (length(CELLLAB_ALL{RXS(rix)}) + 9 + 17):(length(CELLLAB_ALL{RXS(rix)}) + 17 + 17)), 2);

          dIn_pre{rix} = MeanVpre - MeanApre;
          dIn_post{rix} = MeanVpost - MeanApost;
              
          Covres_pre{rix} = cov(residualpre{RXS(rix)}');
          Covres_post{rix} = cov(residualpost{RXS(rix)}');
          Covresinv_pre{rix} = inv(Covres_pre{rix});
          Covresinv_post{rix} = inv(Covres_post{rix});
    
          InputInformation_pre(rix) = dIn_pre{rix}' * Covresinv_pre{rix} * dIn_pre{rix};
          InputInformation_post(rix) = dIn_post{rix}' * Covresinv_post{rix} * dIn_post{rix};
    
    Apre{rix} = xpre{RXS(rix)}(:, 1:size(xpre{RXS(rix)},1));
    Apost{rix} = xpost{RXS(rix)}(:, 1:size(xpost{RXS(rix)},1));
    
    [Vpre{rix}, Dpre{rix}] = eig(Apre{rix});
    [Vpost{rix}, Dpost{rix}] = eig(Apost{rix});
    
    Vinv_pre{rix} = inv(Vpre{rix});
    Vinv_post{rix} = inv(Vpost{rix});
    
    SNRfrac_pre{rix}  = abs(Vinv_pre{rix}  * dIn_pre{rix})  ./ sqrt(diag(Vinv_pre{rix}  * Covres_pre{rix}  * Vinv_pre{rix}'))  / sqrt(InputInformation_pre(rix));
    SNRfrac_post{rix} = abs(Vinv_post{rix} * dIn_post{rix}) ./ sqrt(diag(Vinv_post{rix} * Covres_post{rix} * Vinv_post{rix}')) / sqrt(InputInformation_post(rix));

   
    
    % subselect real eigenmodes
    
    I = find(imag(diag(Dpre{rix})) == 0);
    D = diag(Dpre{rix});
    RealEig_pre{rix} = real(D(I));
    SNRfrac_pre{rix} = real(SNRfrac_pre{rix}(I));
    
    I = find(imag(diag(Dpost{rix})) == 0);
    D = diag(Dpost{rix});
    RealEig_post{rix} = real(D(I));
    SNRfrac_post{rix} = real(SNRfrac_post{rix}(I));

    
end

% combine animals

RealEig_pre_tot = vertcat(RealEig_pre{:});
SNRfrac_pre_tot = vertcat(SNRfrac_pre{:});
RealEig_post_tot = vertcat(RealEig_post{:});
SNRfrac_post_tot = vertcat(SNRfrac_post{:});

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
   
 
    SNRfracbin_pre(i)  = mean(SNRfrac_pre_tot(Ipre( eigs_pre )));
    SNRfracbin_post(i) = mean(SNRfrac_post_tot(Ipost( eigs_post )));

%     SNRfracSEMbin_pre(i)  = std(SNRfrac_pre_tot(Ipre( eigs_pre ))) ./ sqrt(Neigs_pre);
%     SNRfracSEMbin_post(i) = std(SNRfrac_post_tot(Ipost( eigs_post )))  ./ sqrt(Neigs_post);
%        
end


    clear taubin_pre taubin_post taubin_pre_sem taubin_post_sem
    
    for i=1:length(SNRbins)
        
       eigs_pre = and( SNRfrac_pre_tot(Ipre) >= SNRbins(i) - SNRfracbin_width,  SNRfrac_pre_tot(Ipre) < SNRbins(i) + SNRfracbin_width);
       eigs_post = and( SNRfrac_post_tot(Ipost) >= SNRbins(i) - SNRfracbin_width,  SNRfrac_post_tot(Ipost) < SNRbins(i) + SNRfracbin_width);
       Neigs_pre = sum(eigs_pre);
       Neigs_post = sum(eigs_post);              

       taubin_pre(i) = mean(taupre_tot(eigs_pre));
       taubin_post(i) = mean(taupost_tot(eigs_post));

       taubin_pre_sem(i) = std(taupre_tot( eigs_pre)) / sqrt(Neigs_pre);
       taubin_post_sem(i) = std(taupost_tot( eigs_post)) / sqrt(Neigs_post);

    end
    

SNRfracbin_pre_shuff(:,Nshuff) = SNRfracbin_pre;
SNRfracbin_post_shuff(:,Nshuff) = SNRfracbin_post;

taubin_pre_shuff(:,Nshuff) = taubin_pre;
taubin_post_shuff(:,Nshuff) = taubin_post;


edges{1} = 0:10:1750;
edges{2} = 0:0.001:0.25;

sigma_t = 100;
sigma_SNRfrac = 0.025;

s1=(exp(-(bsxfun(@minus, edges{1}, taupre_tot)).^2/(2 * sigma_t^2)));
s2=(exp(-(bsxfun(@minus, edges{2}, SNRfrac_pre_tot(Ipre)).^2)/(2*sigma_SNRfrac^2)));
s3=s1' * s2;
p_pre_smxy_shuff(:,:,Nshuff) = s3;

s1=(exp(-(bsxfun(@minus, edges{1}, taupost_tot)).^2/(2 * sigma_t^2)));
s2=(exp(-(bsxfun(@minus, edges{2}, SNRfrac_post_tot(Ipost)).^2)/(2*sigma_SNRfrac^2)));
s3=s1' * s2;
p_post_smxy_shuff(:,:,Nshuff) = s3;



clearvars -except SNRfracbin_pre_shuff SNRfracbin_post_shuff taubin_pre_shuff taubin_post_shuff p_pre_smxy_shuff p_post_smxy_shuff 

end

edges{1} = 0:10:1750;
edges{2} = 0:0.001:0.25;

for i=1:length(edges{1})
for j=1:length(edges{2})
p_95prc(i,j) = prctile(p_post_smxy_shuff(i,j,:) - p_pre_smxy_shuff(i,j,:),97.5);
p_5prc(i,j)  = prctile(p_post_smxy_shuff(i,j,:) - p_pre_smxy_shuff(i,j,:),2.5);
end
end

for i=1:size(SNRfracbin_pre_shuff,1)
SNRfracbin_post_shuff_5prc(i)  = prctile(SNRfracbin_post_shuff(i,:), 2.5);
SNRfracbin_post_shuff_95prc(i) = prctile(SNRfracbin_post_shuff(i,:), 97.5);
SNRfracbin_pre_shuff_95prc(i)  = prctile(SNRfracbin_pre_shuff(i,:), 97.5); 
SNRfracbin_pre_shuff_5prc(i)   = prctile(SNRfracbin_pre_shuff(i,:), 2.5);

deltaSNRfracbin_shuff_5prc(i)  = prctile(SNRfracbin_post_shuff(i,:) - SNRfracbin_pre_shuff(i,:), 2.5);
deltaSNRfracbin_shuff_95prc(i) = prctile(SNRfracbin_post_shuff(i,:) - SNRfracbin_pre_shuff(i,:), 97.5);
end

for i=1:size(taubin_pre_shuff,1)
taubin_post_shuff_5prc(i) = prctile(taubin_post_shuff(i,:), 2.5);
taubin_post_shuff_95prc(i) = prctile(taubin_post_shuff(i,:), 97.5);
taubin_pre_shuff_95prc(i) = prctile(taubin_pre_shuff(i,:), 97.5);
taubin_pre_shuff_5prc(i) = prctile(taubin_pre_shuff(i,:), 2.5);

deltataubin_shuff_5prc(i) = prctile(taubin_post_shuff(i,:) - taubin_pre_shuff(i,:), 2.5);
deltataubin_shuff_95prc(i) = prctile(taubin_post_shuff(i,:) - taubin_pre_shuff(i,:), 97.5);
end




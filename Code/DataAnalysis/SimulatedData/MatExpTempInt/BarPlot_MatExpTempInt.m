
clear all;

ROOT = '/Users/anguschadwick/Documents/GitHub/SSN_EI_Ring_Models/Code/DataAnalysis/SimulatedData/MatExpTempInt/';
load([ROOT,'TempInt_AllFree.mat'])

Nt = 41

TempIntPre_AllFree = mean(Xpre./Xpre(1,:),1)  * Nt * 0.125;
TempIntPost_AllFree = mean(Xpost./Xpost(1,:),1)  * Nt * 0.125;

[h,p] = ttest(TempIntPre_AllFree,TempIntPost_AllFree, 'tail', 'left')

MagIntPre_AllFree = mean(Xpre) * Nt * 0.125;
MagIntPost_AllFree = mean(Xpost) * Nt * 0.125;


subplot(1,3,1)
hold on
plot([0:0.125:((Nt-1)*0.125)], Xpre(:,1) ./ Xpre(1,1) , 'linewidth', 4, 'color', 'b')
plot([0:0.125:((Nt-1)*0.125)], Xpost(:,1) ./ Xpost(1,1) , 'linewidth', 4, 'color', 'r')
set(gca,'fontsize', 18)
xlabel('Time (s)')
ylabel('Response (dF/F)')
legend('Pre', 'Post')
axis([0,5,0,1])
title('Mouse 1')

subplot(1,3,2)
hold on
shadedErrorBar([0:0.125:((Nt-1)*0.125)], mean(Xpre ./ Xpre(1,:), 2 ), std(Xpre ./ Xpre(1,:),[], 2 ) / sqrt(8), 'lineProps', 'b')
shadedErrorBar([0:0.125:((Nt-1)*0.125)], mean(Xpost ./ Xpost(1,:), 2 ), std(Xpost ./ Xpost(1,:),[], 2 ) / sqrt(8), 'lineProps', 'r')
set(gca,'fontsize', 18)
xlabel('Time (s)')
ylabel('Response (dF/F)')
axis([0,5,0,1])
title('Average over mice')

load([ROOT,'TempInt_FixWeights.mat'])

TempIntPre_FixWeights = TempIntPre  * Nt * 0.125;
TempIntPost_FixWeights = TempIntPost * Nt * 0.125;

load([ROOT,'TempInt_FixInputs.mat'])

TempIntPre_FixInputs = TempIntPre  * Nt * 0.125;
TempIntPost_FixInputs = TempIntPost  * Nt * 0.125;

load([ROOT,'TempInt_AllFixed.mat'])

TempIntPre_AllFixed = TempIntPre  * Nt * 0.125;
TempIntPost_AllFixed = TempIntPost  * Nt * 0.125;

load([ROOT,'TempInt_AplusV.mat'])

TempIntPre_AplusV = TempIntPre  * Nt * 0.125;
TempIntPost_AplusV = TempIntPost  * Nt * 0.125;

load([ROOT,'TempInt_Switching.mat'])

TempIntPre_Switching = TempIntPre  * Nt * 0.125;
TempIntPost_Switching = TempIntPost  * Nt * 0.125;


DataPre_Means = [mean(TempIntPre_AllFree), mean(TempIntPre_FixWeights), mean(TempIntPre_FixInputs), mean(TempIntPre_AllFixed), mean(TempIntPre_AplusV), mean(TempIntPre_Switching)];
DataPost_Means = [mean(TempIntPost_AllFree), mean(TempIntPost_FixWeights), mean(TempIntPost_FixInputs), mean(TempIntPost_AllFixed), mean(TempIntPost_AplusV), mean(TempIntPost_Switching)];

DataPre_SEMs = [std(TempIntPre_AllFree), std(TempIntPre_FixWeights), std(TempIntPre_FixInputs), std(TempIntPre_AllFixed), std(TempIntPre_AplusV), std(TempIntPre_Switching)] / sqrt(8);
DataPost_SEMs = [std(TempIntPost_AllFree), std(TempIntPost_FixWeights), std(TempIntPost_FixInputs), std(TempIntPost_AllFixed), std(TempIntPost_AplusV), std(TempIntPost_Switching)] / sqrt(8);



subplot(1,3,3)
% errorbar([DataPre_Means; DataPost_Means]', [DataPre_SEMs; DataPost_SEMs]')
% axis([0.5, 4.5, 0.2, 0.6])
% set(gca,'fontsize', 18)
% legend('Pre', 'Post')
% ylabel('Time Constant (s)')
% xlabel('Condition')


Means = [DataPre_Means; DataPost_Means]';
SEMs = [DataPre_SEMs; DataPost_SEMs]';
hold on

hb = bar(Means);

ctr =  [0.8571    1.8571    2.8571    3.8571   4.8571  5.8571;    1.1429    2.1429    3.1429    4.1429    5.1429   6.1429]'

hold on
errorbar(ctr', Means', SEMs', '.k', 'linewidth', 3)
hold off
set(gca, 'fontsize', 18)
ylabel('Time Constant (s)')
set(gca, 'xtick', 1:6)
set(gca, 'xticklabel', {'All Free', 'Weights Fixed', 'Inputs Fixed', 'All Fixed', 'Irrelevant Inputs', 'Switching'})
rotateXLabels( gca, 45)


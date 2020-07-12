% close all 
% clear all
% clc
% %% Create or Load Simulated Stimulus
% load('Run4_Data.mat')
NTrain=round(length(tt)*1); % number of training samples
NTest=round(length(tt)*0);% number of test samples
Trial_num=TotCells;
LenWind=1;
%% Test and Train data 
TrainInt=1:NTrain/LenWind;
% TesInt=NTrain/LenWind+1:NTrain/LenWind+NTest/LenWind;
TesInt=TrainInt;
yEstimated=zeros(length(TesInt),1);
TesIntTime=TesInt.*(tt(2)-tt(1));
%% Calculating firing rate and preprocessing the data

 
EE=reshape(Entropy,[45,1]);
 for pa=1:length(WordLengths)
     pa
     for qa=1:length(DeltaT)
         for ps=1:length(WordLengths)
             for qs=1:length(DeltaT)
%        fprintf('p=%d,q=%d \n',p,q);  
    
EEa=resample(EntropyNotNorm{qa,pa}, (DeltaT(qa)/pdt)*WordLengths(pa)*10,10);
EEs=resample(EntropyNotNorm{qs,ps}, (DeltaT(qs)/pdt)*WordLengths(ps)*10,10);
EEs=(EEs-min(EEs))/(max(EEs)-min(EEs));
EEa=(EEa-min(EEa))/(max(EEa)-min(EEa));

fun = @(x)RegressionCost(x,EEs,EEa,normSignal_Mix);
options = optimset('Display','off','MaxIter',100 );

bestx = fminsearch(fun,[0,0],options);
wa=bestx(1);
ws=bestx(2);
EEopt=(wa*EEa+ws*EEs);
rmse(pa,qa,ps,qs)=sqrt(mean(((wa*EEa+ws*EEs)-normSignal_Mix).^2));
                end
             end
     end
 end

%%

[paO,qaO,psO,qsO] = ind2sub(size(rmse), find(rmse == min(rmse(:))));

EEa=resample(EntropyNotNorm{qaO(1),paO(1)}, (DeltaT(qaO(1))/pdt)*WordLengths(paO(1))*10,10);
EEs=resample(EntropyNotNorm{qsO(1),psO(1)}, (DeltaT(qsO(1))/pdt)*WordLengths(psO(1))*10,10);
EEs=(EEs-min(EEs))/(max(EEs)-min(EEs));
EEa=(EEa-min(EEa))/(max(EEa)-min(EEa));

fun = @(x)RegressionCost(x,EEs,EEa,normSignal_Mix);
options = optimset('Display','off','MaxIter',100 );

bestx = fminsearch(fun,[0,0],options);
wa=bestx(1);
ws=bestx(2);
EEopt=(wa*EEa+ws*EEs);
figure
hold on
plot(tt(strInd:endInd),normSignal_Mix(startInd:endInd),'linewidth',2,'color',[0,0,0]);
plot(tt(strInd:endInd),EEopt(startInd:endInd),':r','linewidth',2);
title('Mixed & reconstructed Stimulus','fontsize',22);
xlabel('time (ms)','fontsize',22);
ylabel('Amp.','fontsize',22);
legend({'Mixed stimulus','Reconstructed stimulus'},'fontsize',16);
set (gca, 'fontsize', 15)
hold off
%%

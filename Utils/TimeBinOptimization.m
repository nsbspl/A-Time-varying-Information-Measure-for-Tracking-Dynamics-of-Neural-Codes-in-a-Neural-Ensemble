close all 
clear all
clc
%% Create or Load Simulated Stimulus
% load('Run4_Data.mat')
LenWind=1; % window Length for finding firing rate
NTrain=round(length(tt)*.8); % number of training samples
NTest=round(length(tt)*0.2);% number of test samples
TotCells=30; % number of total neurons
Trial_num=TotCells;
dt=tt(2)-tt(1);
%% Test and Train data 
TrainInt=1:NTrain/LenWind;
TesInt=NTrain/LenWind+1:NTrain/LenWind+NTest/LenWind;
TesIntTime=TesInt.*(tt(2)-tt(1));
%% Calculating firing rate and preprocessing the data
[ySum,yMean,yCells]=MeanFiringRate(F_binary, LenWind); % find mean of firing rate total

[ySumSync,yMeanSync,yCellsSync]=MeanFiringRate(F_binarySync, LenWind); % find mean of firing rate Sync

[ySumAsync,yMeanAsync,yCellsAsync]=MeanFiringRate(F_binaryAsync, LenWind); % find mean of firing rate Asunc

 
 Signal_Mix=conv(Signal_Mix,ones(LenWind,1)/LenWind,'same'); % same length 
 Signal_Mix=downsample(Signal_Mix,LenWind);

 
 Signal_S=conv(Signal_slow(1:end-1)+midpoint,ones(LenWind,1)/LenWind,'same');
 Signal_S=downsample(Signal_S,LenWind);

 
 Signal_F=conv(Signal_fast(1:end-1),ones(LenWind,1)/LenWind,'same');
 Signal_F=downsample(Signal_F,LenWind);
%% Find optimum bin size efor sync spikes

[optNSync, optD, cSync, NSync] = sshist(ySumSync); subplot(3,1,1);plot(NSync,cSync)
fprintf(' Best window length for Sync Spikes is:%d\n',optNSync);
%% Find optimum bin size efor Async spikes
[optNASync, optD, cASync, NASync]= sshist(ySumAsync);  subplot(3,1,2);plot(NASync,cASync)
fprintf(' Best window length for Async Spikes is:%d\n',optNASync);

%% Find optimum bin size efor Mixes spikes
[optNMix, optD, cMix, NMix] = sshist(ySum);  subplot(3,1,3);plot(NMix,cMix)
fprintf(' Best window length for Mixed Spikes is:%d\n',optNMix);
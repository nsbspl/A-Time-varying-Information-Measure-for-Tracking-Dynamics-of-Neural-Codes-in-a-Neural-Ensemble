close all
clear all
clc
%% Create Simulated Stimulus
 run('General_Input.m');
run('Ensemble_Neurons.m') 
% load('Run2_Data.mat')
%% Set  Parameters
 LenWind=100; % window Length for finding firing rate
Cells=1:3; % set of coupled neurons
NTrain=round(length(tt)*0.8); % number of training samples
NTest=round(length(tt)*0.2);% number of test samples
Trial_num=30; % number of total neurons
NumFilters = 10;  % total neumber of regressors 
nthist=20; % spike history length
mode=0;
% a_noise=10; %noise ampl
%% Test and Train data 
TrainInt=1:NTrain/LenWind;
TesInt=int32( NTrain/LenWind+1:NTrain/LenWind+NTest/LenWind);
TesIntTime=int32((TesInt(1)-1)*LenWind+2:length(tt));
%% Calculating firing rate
 [ySum,yMean,yCells]=MeanFiringRate(F_binary, LenWind); % find mean of firing rate total
 Signal_slow=conv(Signal_slow,ones(LenWind,1)/LenWind,'same');
 Signal_slow=downsample(Signal_slow,LenWind);
 %% create the  design matrix
 Xdsgn = DesignMatrix(Signal_Mix,yMean, yCells,mode,NumFilters,nthist,Cells); % create design matrix
%% Coupled and Spike History GLM Possion ASync (Identity)
 yEstimated = GLM_PL(Xdsgn(TrainInt,:), yMean(TrainInt,:),Xdsgn(TesInt,:));
 
%% Result Visualization
yEstimated=resample(yEstimated, LenWind*10,10);
yEstimated(find(yEstimated <0))=0;
yMean= resample(yMean, LenWind*10,10);
figure
subplot(2,1,1)
plot(tt(TesIntTime),yMean(TesIntTime),'linewidth',1);
ylim([0,max(yMean(TesIntTime))]);
title('True Firing Rate');
 xlabel('Time (ms)');
subplot(2,1,2);
plot(tt(TesIntTime), yEstimated,'linewidth',1);
 ylim([0,max(yMean(TesIntTime))]);
title('Estimaed Firing Rate');
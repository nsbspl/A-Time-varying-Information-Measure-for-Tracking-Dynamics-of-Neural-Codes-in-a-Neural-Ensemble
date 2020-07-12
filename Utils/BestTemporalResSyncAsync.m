close all 
clear all
clc
%% Create Simulated Stimulus
 run('General_Input.m');
%% Reduce 1 sample from the end to make the length round 
AsyncStim=Signal_slow(1:end-1,:)+midpoint;
SyncStim=Signal_fast(1:end-1,:);
Signal_Mix=AsyncStim+SyncStim;
tt=tt(1:end-1);
%% parameters definition
TotCells=30; % number of total neurons
Trial_num=TotCells;
a_noise=10; %noise ampl for neurons response
LenWind=1; % window length definition
pdt=(tt(2)-tt(1)); % calculate dt
%% Ensemble_Neurons for ASync stimulus (build async spiking activity)
[F_binaryAsync ]= RM_Ensemble_Neurons(Signal_slow,a_noise, Trial_num,p);
F_binaryAsync=F_binaryAsync(1:end-1,:);
%% Ensemble_Neurons for Sync stimulus (build sync spiking activity)
[F_binarySync ]= RM_Ensemble_Neurons(Signal_fast,a_noise, Trial_num,p);
F_binarySync=F_binarySync(1:end-1,:);
%% Calculate Temporal Resolution for differnt spike types
[ySumSync,yMeanSync,yCellsSync]=MeanFiringRate(F_binarySync, LenWind); % find mean of firing rate Sync
[ySumAsync,yMeanAsync,yCellsAsync]=MeanFiringRate(F_binaryAsync, LenWind); % find mean of firing rate Asunc
%%
clc
bestSigmaSync=FindBestTemporalResolution(ySumSync, false)
bestSigmaASync=FindBestTemporalResolution(ySumAsync, false)
% bestSigmaMixed=FindBestTemporalResolution(ySum, false)
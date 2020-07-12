close all
clear all
clc
addpath('Utils')
%% Create Data 
% parameter setting
SamplingTime = .05;% msec
EndTime = 10;%sec
TotCells=100; % number of total neurons
Trial_num=TotCells;
a_noise=30; %noise ampl

fprintf('Generating Dataset  ... \n');
run('CreateDataSet.m');

fprintf('Dataset finished,SamplingTime=%2.2f ms,Length= %3.3f s,NumberofCells=%d,noiseAmp=%2.2f \n', ...
   SamplingTime, EndTime, TotCells, a_noise );

%% Deta Visualization
fprintf('Running Data Visualization ... \n');
run('CreateFigers.m');
fprintf(' is done\n');
fprintf('See figures 1- 5 \n');
%% Run distribution estimation
fprintf('Running distribution estimation ... \n');
run('Distribution_estimation_nonparametric.m');
fprintf('is done\n');
fprintf('See figure 6  \n');
%% Run direct information estimation from spiking activity
fprintf('Running direct information estimation from spiking activity ... \n');
run('EntropyAsFunctionOfTemporalAndStructuralRes.m');
fprintf(' is done\n');
fprintf('See figures 7 & 8 \n');
%% Run Kalman-Filter decoder
fprintf('Running KF-Decoder ... \n');
run('KFDecoderModel.m');
fprintf(' is done\n');
fprintf('See figures 9-11\n');

%% Run TVE measure
fprintf('Running TVE Measure ... \n');
run('TvEntropySpikeIdentification.m');
fprintf('TVE Measure is done\n');
fprintf('See rest of figures\n');
%% simple regression
fprintf('find best setting for TVE to reconstruct the stimulation ... \n');
run('Decoder_Entrpy.m');
fprintf(' is done\n');
fprintf('See the last figure\n');

%% Similarty between TVI and spectrum
fprintf('Similarty between TVI and spectrum ... \n');
run('TVE_Spec.m');
fprintf(' is done\n');
fprintf('See the last figure\n');


%%%%%%%%% The End %%%%%%%%%%%%%%%%%%%%%

%% Create Simulated Stimulus
addpath('IEEE_CIM')
addpath('Info_Transmission')
addpath('Information Theory and Neuroscience')

% I sweeped two next line of codes
p = NeuronType(SamplingTime,EndTime,'CD');
tt = 0:p.dt:EndTime*1e3;
p.V3 = -15;%-19;% + randn(1);%Integ_c(1,end);
p.exc_avg = 0.8;%1.12;%Integ_c(2,end);
p.inh_avg = 1.7;

stimDist = 'Gauss';
tEnd = 0.05; % Total synaptic duation (sec)
tau_rise = 0.5; % msec
tau_fall = 3; % msec --> stand for AMPA (Exc) synapses
Gz = exp2syn_new(p.dt,tEnd,tau_rise,tau_fall);
L = EndTime/p.dt*1e+3 + 1;
Sigma = [0.1; 0.2; 0.3; 0.4; 0.5];
amp_thresh = 0.031;
run('General_Input.m');

%% Reduce 1 sample from end
AsyncStim=Signal_slow(1:end-1,:)+midpoint;
SyncStim=Signal_fast(1:end-1,:);
Signal_Mix=AsyncStim+SyncStim;
tt=tt(1:end-1);
%% Set  Parameters

%% Ensemble_Neurons for ASync stimulus
[F_binaryAsync ]= RM_Ensemble_Neurons(Signal_slow+midpoint,a_noise, Trial_num,p);
F_binaryAsync=F_binaryAsync(1:end-1,:);
%% Ensemble_Neurons for Sync stimulus
[F_binarySync ]= RM_Ensemble_Neurons(Signal_fast,a_noise, Trial_num,p);
F_binarySync=F_binarySync(1:end-1,:);
%% Ensemble_Neurons for mixed stimulus
[F_binary ]= RM_Ensemble_Neurons(Signal_slow+midpoint+Signal_fast,a_noise, Trial_num,p);
F_binary=F_binary(1:end-1,:);
%% Ensemble_Neurons for no stimulus
F_binary_noise=RM_Ensemble_Neurons(Signal_slow*0,a_noise, Trial_num,p);
F_binary_noise=F_binary_noise(1:end-1,:);
Signal_fast=Signal_fast(1:end-1,:);
Signal_slow=Signal_slow(1:end-1,:);
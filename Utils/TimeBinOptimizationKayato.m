
%% Create or Load Simulated Stimulus

dt=tt(2)-tt(1);
LenWind=1;
TotCells=1;
%% Calculating firing rate and preprocessing the data
[ySum,yMean,yCells]=MeanFiringRate(F_binary, LenWind); % find mean of firing rate total

[ySumSync,yMeanSync,yCellsSync]=MeanFiringRate(F_binarySync, LenWind); % find mean of firing rate Sync

[ySumAsync,yMeanAsync,yCellsAsync]=MeanFiringRate(F_binaryAsync, LenWind); % find mean of firing rate Asunc

%% Find optimum bin size efor sync spikes
N=1:1:200;  
[optNSync,cSync, NSync] = OptimizeTimeBin(ySumSync,N,TotCells,dt); subplot(3,1,1);plot(NSync,cSync)
fprintf(' Best dt for Sync Spikes is:%d\n',optNSync*dt);
%% Find optimum bin size efor Async spikes
% N=1:1:10; 
[optNASync,cASync, NASync] = OptimizeTimeBin(ySumAsync,N,TotCells,dt);  subplot(3,1,2);plot(NASync,cASync)
fprintf(' Best dt for Async Spikes is:%d\n',optNASync*dt);

%% Find optimum bin size efor Mixes spikes
% N=1:1:10; 
[optNMix,cMix, NMix] = OptimizeTimeBin(ySum,N,TotCells,dt);  subplot(3,1,3);plot(NMix,cMix)
fprintf(' Best dt for Mixed Spikes is:%d\n',optNMix*dt);
% close all 
% clear all
% clc
%% Create Simulated Stimulus
% load('Run4_Data.mat')
% load('Run4_Data_noise.mat')
% F_binary_noise=F_binary_noise(1:end-1,:);
Max=10;
%%
WordLengths=[25, 20,10, 5, 4, 2,1 ]; % word length (structure resolutions)
DeltaTs=[.05, .1, .2, .5, 1, 10, 20,100]; % temporal resolutions
pdt=(tt(2)-tt(1)); % dt
%% Calculate Entropy as function of  
% window with length L (structure resolution
% and temporal resolution) for all spikes
EntropyTotal = entropyFuncTempResolutionAndWindoL(F_binary_noise(:,1:Max), DeltaTs,WordLengths,pdt);

%% Calculate Entropy as function of  
% window with length L (structure resolution
% and temporal resolution) for Async spikes
EntropyAsync= entropyFuncTempResolutionAndWindoL(F_binaryAsync(:,1:Max), DeltaTs,WordLengths,pdt);

%% Calculate Entropy as function of  
% window with length L (structure resolution
% and temporal resolution) for Sync spikes
EntropySync = entropyFuncTempResolutionAndWindoL(F_binarySync(:,1:Max), DeltaTs,WordLengths,pdt);

%% Calculate Entropy as function of  
% window with length L (structure resolution
% and temporal resolution) for Sync spikes
EntropyMixed = entropyFuncTempResolutionAndWindoL(F_binary(:,1:Max), DeltaTs,WordLengths,pdt);
%% plot
figure
subplot(2,2,1)
clims = ([min(min(EntropyMixed)),max(max(EntropyMixed))]);
imagesc(1./WordLengths,log(DeltaTs),(EntropyTotal),clims)
title('Entropy of no-stimulus spikes','fontsize',22);
xlabel('1/L','fontsize',22);
ylabel('log(dt)','fontsize',22);
set (gca, 'fontsize', 15)
% colorbar
subplot(2,2,2)

imagesc(1./WordLengths,log(DeltaTs),(EntropyAsync),clims)
title('Entropy of Async-Spikes','fontsize',22);
xlabel('1/L','fontsize',22);
ylabel('log(dt)','fontsize',22);
set (gca, 'fontsize', 15)
% zlim([min(min(EntropyMixed)),max(max(EntropyMixed))]);
% colorbar
subplot(2,2,3)

imagesc(1./WordLengths,log(DeltaTs),(EntropySync),clims)
title('Entropy of Sync-Spikes','fontsize',22);
xlabel('1/L','fontsize',22);
ylabel('log(dt)','fontsize',22);
set (gca, 'fontsize', 15)
% zlim([min(min(EntropyMixed)),max(max(EntropyMixed))]);
% colorbar
subplot(2,2,4)

imagesc(1./WordLengths,log(DeltaTs),(EntropyMixed),clims)
title('Entropy of Mixed-Spikes','fontsize',22);
xlabel('1/L','fontsize',22);
ylabel('log(dt)','fontsize',22);
set (gca, 'fontsize', 15)
% zlim([min(min(EntropyMixed)),max(max(EntropyMixed))]);
figure
imagesc(1./WordLengths,log(DeltaTs),(EntropyMixed),clims)
colorbar()
%%
figure
subplot(411)
indxs=1;
plot(1./WordLengths,(EntropyTotal(7,:)),'linewidth',2)
title('Entropy no-stimulus','fontsize',22);
% xlabel('1/L');
% ylabel('En. bit/sec');
set (gca, 'fontsize', 15)
subplot(412)
plot(1./WordLengths,(EntropyAsync(7,:)),'linewidth',2)
title('Entropy Async Stimulus','fontsize',22);
% xlabel('1/L');
% ylabel('En');
set (gca, 'fontsize', 15)
subplot(413)
plot(1./WordLengths,(EntropySync(7,:)),'linewidth',2)
title('Entropy Sync Stimulus','fontsize',22);
% xlabel('1/L');
% ylabel('En');
set (gca, 'fontsize', 15)

subplot(414)
plot(1./WordLengths,(EntropyMixed(7,:)),'linewidth',2)
title('Entropy Mixed Stimulus','fontsize',22);
xlabel('1/L','fontsize',22);
ylabel('En. bit/sec','fontsize',22);
% legend(['dt=', num2str(DeltaTs(indxs(1)))],['dt=', num2str(DeltaTs(indxs(2)))],['dt=', num2str(DeltaTs(indxs(3)))])
set (gca, 'fontsize', 15)
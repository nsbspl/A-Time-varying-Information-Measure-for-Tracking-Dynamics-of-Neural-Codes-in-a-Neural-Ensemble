% close all 
% clear all
% clc

%% Calculating firing rate and preprocessing the data
[ySum,yMean,yCells]=MeanFiringRate(F_binary, LenWind); % find mean of firing rate total
yMean=(yMean-min(yMean))/(max(yMean)-min(yMean));

[ySumSync,yMeanSync,yCellsSync]=MeanFiringRate(F_binarySync, LenWind); % find mean of firing rate Sync
yMeanSync=(yMeanSync-min(yMean))/(max(yMean)-min(yMean));

[ySumAsync,yMeanAsync,yCellsAsync]=MeanFiringRate(F_binaryAsync, LenWind); % find mean of firing rate Asunc
 yMeanAsync=(yMeanAsync-min(yMean))/(max(yMean)-min(yMean));
 
%% Distribution Estimation of spikes
bandWidth=.4
figure
%%%%%%%%%%%%% Sync %%%%%%%%%
subplot(211)
syncDist='Binomial'; % distribution type
xSync=min(ySum):1:max(ySum);
hSync=hist(ySumSync,length(xSync)); % actual distribution
hSync=hSync./sum(hSync);
dPSync=fitdist(ySumSync,'Kernel','BandWidth',bandWidth ); % parametric estimation of the distribution
ySync = pdf(dPSync,xSync);
hold on % compare the results
bar(xSync,hSync,0.5,'b');
bar(xSync,ySync,0.25,'r');
ylabel('Density','fontsize',22);
xlim([min(ySum),max(ySum)]);
legend({'true','estimated'},'fontsize',16);
hold off
R2=r_squre(ySync(1:end),hSync(1:end));
%title(['Distribution Sync Spikes, dist=',syncDist,' with lambda=', num2str(dPSync.lambda)],'fontsize',22);
set (gca, 'fontsize', 15)
%%%%%%%%%%%%% Async %%%%%%%%%
subplot(212)
asyncDist='poisson';
xAsync=min(ySum):1:max(ySum);
hAsync=hist(ySumAsync,length(xAsync));
hAsync=hAsync./sum(hAsync);
dPAsync=fitdist(ySumAsync,'Kernel','BandWidth',bandWidth);
yAsync = pdf(dPAsync,xAsync);
hold on
bar(xAsync,hAsync,0.5,'b');
bar(xAsync,yAsync,0.25,'r');
legend({'true','estimated'},'fontsize',16);
xlim([min(ySum),max(ySum)]);
hold off
R2=r_squre(yAsync,hAsync);
%title(['Distribution Async Spikes, dist=',asyncDist,' with lambda=', num2str(dPAsync.lambda)],'fontsize',22);
set (gca, 'fontsize', 15)

%%
fprintf('KL-Divergence Between   Sync  & Async spikes distributions =%f \n',KLDiv(hSync,hAsync));
% fprintf('KL-Divergence Between distribution and Poisson estimated  of Sync =%f \n',KLDiv(ySync,hSync));
% fprintf('KL-Divergence Between distribution and Poisson estimated  of Async =%f \n',KLDiv(yAsync,hAsync));
fprintf('KL-Divergence Between Poisson Estimation of  Sync  & Async spikes distributions =%f \n',KLDiv(ySync,yAsync));
% fprintf('KL-Divergence Between distribution and normal estimated  of Async =%f \n',KLDiv(yNAsync,hAsync));
% fprintf('KL-Divergence Between distribution and binomial estimated  of Sync =%f \n',KLDiv(yBSync,hSync));
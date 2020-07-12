% close all
% clear all
% clc
%%
%load('R5.mat')
%%

pp=2;
windowLength=150;
[ySum,yMean,yCells]=MeanFiringRate(F_binary, windowLength);
yMean=resample(yMean, windowLength*10,10);
ySum=resample(ySum, windowLength*10,10);
Spect=zeros(length(DeltaT),length(startInd:endInd));
SpectF=zeros(length(DeltaT),length(startInd:endInd));
SpectS=zeros(length(DeltaT),length(startInd:endInd));
 for q=1:length(DeltaT)
%        fprintf('p=%d,q=%d \n',p,q); 
% En=(EntropyNotNorm{q,pp}-min(EntropyNotNorm{q,pp}))/(max(EntropyNotNorm{q,pp})-min(EntropyNotNorm{q,pp}));
% temp=resample(EntropyNotNorm{q,pp}, (DeltaT(q)/pdt)*WordLengths(pp)*10,10);
temp=Entropy{q,pp};

 if q <=5 
     SpectF(q,:)=(abs(temp(startInd:endInd)));
 else
     SpectS(q,:)=(abs(temp(startInd:endInd)));
 end
 Spect(q,:)=SpectS(q,:)+SpectF(q,:);
 end
gausWindowLength=windowLength;
gausWindowSigma=windowLength/1;
gausWindow = fspecial('gaussian', gausWindowLength, gausWindowSigma);
gausWindow=mean(gausWindow);
gausWindow=gausWindow/max(gausWindow);
[ySums,yMeans, yCellss ]= smoothNeuralActivity_v2(Spect', gausWindow);
%% F1
figure
subplot(4,1,1:2)
imagesc(tt(strInd:endInd),log10(DeltaT(1:end)), (yCellss(:,1:end))');
title('TVE Spectrum','fontsize',22);
set (gca, 'fontsize', 15)
ylabel('log(dt)','fontsize',22);
subplot(413)
plot(tt(startInd:endInd),yMean(startInd:endInd),'linewidth',2,'color',[0,0,0]);
title('Mixed firing rate','fontsize',22);
ylim([0,1]);
set (gca, 'fontsize', 15)
subplot(414)
plot(tt(startInd:endInd),normSignal_Mix(startInd:endInd),'linewidth',2,'color',[0,0,0]);
title('Mixed Stimulus','fontsize',22);
xlabel('time (ms)','fontsize',22);
ylabel('Amp.','fontsize',22);
set (gca, 'fontsize', 15)
%% F2

H=fspecial('disk',5);
yCellss = imfilter(yCellss,H,'replicate');

figure
subplot(4,1,1:2)
imagesc(tt(strInd:endInd),log10(DeltaT(1:end)), (yCellss(:,1:end))');
title('TVE Spectrum','fontsize',22);
set (gca, 'fontsize', 15)
ylabel('log(\delta{t})','fontsize',22);
subplot(413)
plot(tt(startInd:endInd),yMean(startInd:endInd),'linewidth',2,'color',[0,0,0]);
title('Mixed firing rate','fontsize',22);
ylim([0,1]);
set (gca, 'fontsize', 15)
subplot(414)
plot(tt(startInd:endInd),normSignal_Mix(startInd:endInd),'linewidth',2,'color',[0,0,0]);
title('Mixed Stimulus','fontsize',22);
xlabel('time (ms)','fontsize',22);
ylabel('Amp.','fontsize',22);
set (gca, 'fontsize', 15)

%% F3
Hs=fspecial('disk',5);
Hf=fspecial('disk',2);
[ySums,yMeans, yCellssF ]= smoothNeuralActivity_v2(SpectF', gausWindow);
yCellssF = imfilter(yCellssF,Hf,'replicate');
[ySums,yMeans, yCellssS ]= smoothNeuralActivity_v2(SpectS', gausWindow);
yCellssS = imfilter(yCellssS,Hs,'replicate');
figure
clims = ([min(yCellss(:)),max(yCellss(:))]);
subplot(311)
imagesc(tt(strInd:endInd),log10(DeltaT(1:end)), (yCellssF(:,1:end))',clims);
title('TVE Spectrum Sync-spikes','fontsize',22);
ylabel('log(\delta{t})','fontsize',22);
set (gca, 'fontsize', 15)
subplot(312)
imagesc(tt(strInd:endInd),log10(DeltaT(1:end)), (yCellssS(:,1:end))',clims);
title('TVE Spectrum Async-spikes','fontsize',18);
ylabel('log(dt)','fontsize',22);
set (gca, 'fontsize', 15)
subplot(313)
imagesc(tt(strInd:endInd),log10(DeltaT(1:end)), (yCellss(:,1:end))',clims);
title('TVE Spectrum Mixed','fontsize',22);
ylabel('log(dt)','fontsize',22);
set (gca, 'fontsize', 15)
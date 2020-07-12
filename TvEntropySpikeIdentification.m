% close all 
% clear all
% clc
%% Create or Load Simulated Stimulus
% load('Run4_Data.mat')
%%
WordLengths=[  20,10,5, 2,1 ]; % word lengths
pdt=SamplingTime; % calculate dt
DeltaT=[ .05, .1, .2, .5, 1, 5, 10, 50, 100]; % temporal resolutions
startInd=1;
endInd=length(F_binary);
%% Temporal resolution

[signalLength, NumCells] = size(F_binary(startInd:endInd,:));
normSignal_Mix=(Signal_Mix(startInd:endInd,:)-min(Signal_Mix(startInd:endInd,:)))/(max(Signal_Mix(startInd:endInd,:))-min(Signal_Mix(startInd:endInd,:)));
normSignal_fast=(Signal_fast(startInd:endInd,:)-min(Signal_fast(startInd:endInd,:)))/(max(Signal_fast(startInd:endInd,:))-min(Signal_fast(startInd:endInd,:)));
normSignal_slow=(Signal_slow(startInd:endInd,:)-min(Signal_slow(startInd:endInd,:)))/(max(Signal_slow(startInd:endInd,:))-min(Signal_slow(startInd:endInd,:)));
corResult=zeros(length(DeltaT),length(WordLengths),2);
q=0;
Entropy={};
EntropyNotNorm={};
for dt = DeltaT
    q=q+1;
    BinLength = (dt/pdt);
    Sp=zeros(signalLength / BinLength,NumCells);
    for i = 1: signalLength / BinLength
        for p=1:NumCells
            if sum( F_binary( (i-1)*BinLength+1: (i)*BinLength,p)) >=1
                Sp(i,p) =1;
            end
        end
    end
    %% Calculate Entropy as function of  
% window with length L 
% (dt is constant here and is the maximium resolution we can use)

p=0;
for i =WordLengths
    p=p+1;
    Entropy{q,p} = entropyFuncTimeIndexWindoL(Sp,i,dt);
end
%% Pre Plot (make the graphs smoother)
p=0;
for i =WordLengths
    p=p+1;
    EntropyNotNorm{q,p}=Entropy{q,p};
    Entropy{q,p} =resample(Entropy{q,p}, BinLength*i*10,10);
    Entropy{q,p}=(Entropy{q,p}-min(Entropy{q,p}))/(max(Entropy{q,p})-min(Entropy{q,p}));
    corResult(q,p,1)=sqrt(mean((normSignal_fast-Entropy{q,p}).^2));
    corResult(q,p,2)=sqrt(mean((normSignal_slow-Entropy{q,p}).^2));
end
%% Plot the result
figure

subplot(length(WordLengths)+1,1,1);
 plot(Signal_Mix);
 ylim([min(Signal_Mix),max(Signal_Mix)]);
 title('signal Mixed','fontsize',22);
 for i = 1 : length(WordLengths)
  subplot(length(WordLengths)+1,1,i+1);
 plot(Entropy{q,i});
%  ylim([0,max(Entropy{i})]);
 title(['Entropy for word Length=', num2str(WordLengths(i)),',timeRes=',num2str(dt)],'fontsize',22);
 end
 xlabel('time (ms)','fontsize',22);
end
set (gca, 'fontsize', 15)
%% plot corrolation signals
startInd=1;
endInd=20*1000*5;
corResult=[];

for p=1:length(WordLengths)
    for q=1:length(DeltaT)
     corResult(q,p,1)=corr(normSignal_fast(startInd:endInd),Entropy{q,p}(startInd:endInd));
    corResult(q,p,2)=corr(normSignal_slow(startInd:endInd),Entropy{q,p}(startInd:endInd));
    corResult(q,p,3)=corr(normSignal_Mix(startInd:endInd),Entropy{q,p}(startInd:endInd));
    end
end
H=fspecial('disk',5);
corResult = imfilter(corResult,H,'replicate');
clims = ([min(min(corResult(:,:,3))),max(max(corResult(:,:,3)))]);
figure
subplot(3,1,1)
imagesc(1./WordLengths, log10(DeltaT), corResult(:,:,1),clims)
xlabel('1/L','fontsize',22);
ylabel( 'log(dt)','fontsize',22);
title('TVE corr. with fast signal','fontsize',22);
set (gca, 'fontsize', 15)
subplot(3,1,2)
imagesc(1./WordLengths, log10(DeltaT),corResult(:,:,2),clims)
title('TVE corr. with slow signal','fontsize',22);
xlabel('1/L','fontsize',22);
ylabel( 'log(dt)','fontsize',22);

set (gca, 'fontsize', 15)
subplot(3,1,3)
imagesc(1./WordLengths, log10(DeltaT),corResult(:,:,3),clims)
title('TVE corr. with Mixed signal','fontsize',22);
xlabel('1/L','fontsize',22);
ylabel( 'log(dt)','fontsize',22);
set (gca, 'fontsize', 15)
figure
imagesc(1./WordLengths, log10(DeltaT),corResult(:,:,3),clims)
colorbar('southoutside')
%% Sync

figure
subplot(4,1,1)
plot((startInd:endInd).*pdt,normSignal_fast(startInd:endInd),'linewidth',2,'color',[1,0,0]);
% xlabel('Time (ms)');
ylabel('Amp.','fontsize',22);
title('Fast Stimulus','fontsize',22);
set (gca, 'fontsize', 15)
subplot(4,1,2)
q=1;
p=5;
plot((startInd:endInd).*pdt,Entropy{q,p}(startInd:endInd),'linewidth',2,'color',[0,0,1]);
% xlabel('Time (ms)');
ylabel('Amp.','fontsize',22);
title(['TVE measure for word Length=', num2str(WordLengths(p)),',timeRes=',num2str(DeltaT(q))],'fontsize',22);
set (gca, 'fontsize', 15)
subplot(4,1,3)
q=1;
p=2;
plot((startInd:endInd).*pdt,Entropy{q,p}(startInd:endInd),'linewidth',2,'color',[0,0,1]);
% xlabel('Time (ms)');
ylabel('Amp.','fontsize',22);
title(['TVE measure for word Length=', num2str(WordLengths(p)),',timeRes=',num2str(DeltaT(q))],'fontsize',22);
set (gca, 'fontsize', 15)
subplot(4,1,4)
q=6;
p=2;
plot((startInd:endInd).*pdt,Entropy{q,p}(startInd:endInd),'linewidth',2,'color',[0,0,1]);
title(['TVE measure for word Length=', num2str(WordLengths(p)),',timeRes=',num2str(DeltaT(q))],'fontsize',22);
xlabel('Time (ms)','fontsize',22);
ylabel('Amp.','fontsize',22);
set (gca, 'fontsize', 15)
%% Async
figure
subplot(4,1,1)
plot((startInd:endInd).*pdt,normSignal_slow(startInd:endInd),'linewidth',2,'color',[1,0,0]);
% xlabel('Time (ms)');
ylabel('Amp.','fontsize',22);
title('Slow Stimulus','fontsize',22);
set (gca, 'fontsize', 15)

subplot(4,1,2)
q=1;
p=5;
plot((startInd:endInd).*pdt,Entropy{q,p}(startInd:endInd),'linewidth',2,'color',[0,0,1]);
% xlabel('Time (ms)');
ylabel('Amp.','fontsize',22);
title(['TVE measure for word Length=', num2str(WordLengths(p)),',timeRes=',num2str(DeltaT(q))],'fontsize',22);
set (gca, 'fontsize', 15)
subplot(4,1,3)
q=7;
p=2;
plot((startInd:endInd).*pdt,Entropy{q,p}(startInd:endInd),'linewidth',2,'color',[0,0,1]);
% xlabel('Time (ms)');
ylabel('Amp.','fontsize',22);
title(['TVE measure for word Length=', num2str(WordLengths(p)),',timeRes=',num2str(DeltaT(q))],'fontsize',22);
set (gca, 'fontsize', 15)
subplot(4,1,4)
q=9;
p=2;
plot((startInd:endInd).*pdt,Entropy{q,p}(startInd:endInd),'linewidth',2,'color',[0,0,1]);
title(['TVE measure for word Length=', num2str(WordLengths(p)),',timeRes=',num2str(DeltaT(q))],'fontsize',22);
xlabel('Time (ms)','fontsize',22);
ylabel('Amp.','fontsize',22);
set (gca, 'fontsize', 15)
%%
%% Sync & Async
figure
subplot(6,1,1)
plot((startInd:endInd).*pdt,normSignal_Mix(startInd:endInd),'linewidth',2,'color',[1,0,0]);
% xlabel('Time (ms)');
% ylabel('Amp.');
title('Mixed Stimulus','fontsize',14);
% ax = gca; % current axes
% ax.FontSize = 12;
set (gca, 'fontsize', 15)
subplot(6,1,2)
q=1;
p=5;
plot((startInd:endInd).*pdt,Entropy{q,p}(startInd:endInd),'linewidth',2,'color',[0,0,1]);
% xlabel('Time (ms)');
% ylabel('Amp.');
title(['TVE measure for word Length=', num2str(WordLengths(p)),',timeRes=',num2str(DeltaT(q))],'fontsize',14);
% ax = gca; % current axes
% ax.FontSize = 12;
set (gca, 'fontsize', 15)
subplot(6,1,3)
qs=1;
ps=2;
plot((startInd:endInd).*pdt,Entropy{qs,ps}(startInd:endInd),'linewidth',2,'color',[0,0,1]);
% xlabel('Time (ms)');
% ylabel('Amp.');
title(['TVE measure for word Length=', num2str(WordLengths(ps)),',timeRes=',num2str(DeltaT(qs))],'fontsize',14);
% ax = gca; % current axes
% ax.FontSize = 12;
set (gca, 'fontsize', 15)
subplot(6,1,4)
qa=6;
pa=2;
plot((startInd:endInd).*pdt,Entropy{qa,pa}(startInd:endInd),'linewidth',2,'color',[0,0,1]);
% xlabel('Time (ms)');
% ylabel('Amp.');
title(['TVE measure for word Length=', num2str(WordLengths(pa)),',timeRes=',num2str(DeltaT(qa))],'fontsize',14);
% ax = gca; % current axes
% ax.FontSize = 12;
set (gca, 'fontsize', 15)
subplot(6,1,5)
q=7;
p=2;
plot((startInd:endInd).*pdt,Entropy{q,p}(startInd:endInd),'linewidth',2,'color',[0,0,1]);
title(['TVE measure for word Length=', num2str(WordLengths(p)),',timeRes=',num2str(DeltaT(q))],'fontsize',14);
% xlabel('Time (ms)','fontsize',14);
% ylabel('Amp.');
% ax = gca; % current axes
% ax.FontSize = 12;
set (gca, 'fontsize', 15)
subplot(6,1,6)
q=9;
p=2;
plot((startInd:endInd).*pdt,Entropy{q,p}(startInd:endInd),'linewidth',2,'color',[0,0,1]);
title(['TVE measure for word Length=', num2str(WordLengths(p)),',timeRes=',num2str(DeltaT(q))],'fontsize',14);
xlabel('time (ms)','fontsize',22);
set (gca, 'fontsize', 15)
%%

%% marginalize
margEntropy=zeros(size(EntropyNotNorm));
[M,N]=size(EntropyNotNorm);

for p=1:M
    for q=1:N
        margEntropy(p,q)=2*mean(EntropyNotNorm{p,q});
    end
end
figure
H=fspecial('disk',10);
margEntropy = imfilter(margEntropy,H,'replicate');

imagesc(1./WordLengths,log(DeltaT),(margEntropy))
title('Mean TVE over time','fontsize',22);
xlabel('1/L','fontsize',22);
ylabel('log(dt)','fontsize',22);
colorbar
set (gca, 'fontsize', 15)
% close all 
% clear all
% clc
%% Create Simulated Stimulus

% load('Run4_Data.mat');

%% Set  Parameters
 LenWind=100; % window Length for finding firing rate

 %% First figure
 figure
 strInd=1;
 endInd=20*1000*5;
 Signal_SS=Signal_slow;
 subplot(2,1,1)
 plot(tt(strInd:endInd) ,Signal_SS(strInd:endInd),'-','linewidth',2,'color',[0,0,1]);
 title('Slow Stimulus','fontsize',22);
 xlabel('time (ms)','fontsize',22);
 ylabel('Amp.','fontsize',22);
 Signal_FS=Signal_fast;
set (gca, 'fontsize', 15)
 subplot(2,1,2)
 plot(tt(strInd:endInd) ,Signal_FS(strInd:endInd),'linewidth',2,'color',[1,0,0]);
 title('Fast Stimulus','fontsize',22);
 xlabel('time (ms)','fontsize',22);
 ylabel('Amp.','fontsize',22);
set (gca, 'fontsize', 15)
 %% Second Figure
 figure
 subplot(2,1,1)
 plot(tt(strInd:endInd),Signal_SS(strInd:endInd),'linewidth',2,'color',[0,0,1]);
 hold on
 FastIndx=find(Signal_FS(strInd:endInd) ~= 0);
 plot(tt(FastIndx) ,Signal_FS(FastIndx)+Signal_SS(FastIndx),'.','linewidth',2,'color',[1,0,0]);
 hold off
 title('Mixing Stimulus','fontsize',22);
 xlabel('time (ms)','fontsize',22);
 ylabel('Amp.','fontsize',22);
 legend({'Slow Stimulus','Fast Stimulus'},'FontSize',16);
set (gca, 'fontsize', 15)
 subplot(2,1,2)
 plot(tt(strInd:endInd) ,Signal_SS(strInd:endInd)+Signal_FS(strInd:endInd),'linewidth',2,'color',[0,0,0]);
 
 title('Mixed Stimulus','fontsize',22);
 xlabel('time (ms)','fontsize',22);
 ylabel('Amp.','fontsize',22);
 set (gca, 'fontsize', 15)
 %% Third Figure
 
 figure
 subplot(2,1,1)
 hold on
 sz = 80;
 for i=1:TotCells
 scatter(tt(strInd:endInd),F_binaryAsync(strInd:endInd,i)*i,sz,'.b');
 end
 hold off
 title('Async. Spikes','fontsize',22);
 xlabel('time (ms)','fontsize',22);
 ylabel('neurons','fontsize',22);
 ylim([.5,TotCells+.5]);
  set (gca, 'fontsize', 15)
  subplot(2,1,2);
 hold on
 for i=1:TotCells
 scatter(tt(strInd:endInd) ,F_binarySync(strInd:endInd,i)*i,sz,'.r');
 end
 hold off
 title('Sync. Spikes','fontsize',22);
 xlabel('time (ms)','fontsize',22);
 ylabel('neurons','fontsize',22);
 ylim([.5,TotCells+.5]);
 set (gca, 'fontsize', 15)
 %% Fourth Figure
  figure
 subplot(2,1,1)
 hold on
 sz = 80;
 for i=1:TotCells
 scatter(tt(strInd:endInd),F_binaryAsync(strInd:endInd,i)*i,sz,'.b');
 scatter(tt(strInd:endInd) ,F_binarySync(strInd:endInd,i)*i,sz,'.r');
 end
 hold off
 title('Mixing Spikes','fontsize',22);
 xlabel('time (ms)','fontsize',22);
 ylabel('neurons','fontsize',22);
 ylim([.5,TotCells+.5])
  legend({'Async. Spikes','Sync. Spikes'},'FontSize',16);
  set (gca, 'fontsize', 15)
 subplot(2,1,2)
 hold on
 for i=1:TotCells
 scatter(tt(strInd:endInd) ,F_binary(strInd:endInd,i)*i,sz,'.k');
 end
 hold off
 title('Mixed Spikes','fontsize',22);
 xlabel('time (ms)','fontsize',22);
 ylabel('neurons','fontsize',22);
 ylim([.5,TotCells+.5])
 set (gca, 'fontsize', 15)
 %%
  figure
 subplot(2,1,1)
 plot(tt(strInd:endInd),Signal_SS(strInd:endInd),'linewidth',2,'color',[0,0,1]);
 hold on
 FastIndx=find(Signal_FS(strInd:endInd) ~= 0);
 plot(tt(FastIndx) ,Signal_FS(FastIndx)+Signal_SS(FastIndx),'.','linewidth',2,'color',[1,0,0]);
 hold off
 title('Mixing Stimulus','fontsize',22);
 xlabel('time (ms)','fontsize',22);
 ylabel('Amp.','fontsize',22);
 legend({'Slow Stimulus','Fast Stimulus'},'FontSize',16)
 set (gca, 'fontsize', 15)
 subplot(2,1,2)
 hold on
 sz = 80;
 for i=1:TotCells
 scatter(tt(strInd:endInd),F_binaryAsync(strInd:endInd,i)*i,sz,'.b');
 scatter(tt(strInd:endInd) ,F_binarySync(strInd:endInd,i)*i,sz,'.r');
 end
 hold off
 title('Mixing Spikes','fontsize',22);
 xlabel('time (ms)','fontsize',22);
 ylabel('neurons','fontsize',22);
 ylim([.5,TotCells+.5]);
  legend({'Async. Spikes','Sync. Spikes'},'FontSize',16);
 set (gca, 'fontsize', 15)
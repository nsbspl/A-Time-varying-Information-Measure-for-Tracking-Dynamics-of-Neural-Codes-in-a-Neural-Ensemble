% close all 
% clear all
% clc
% %% Create or Load Simulated Stimulus
% load('Run4_Data.mat')
LenWind=100; % window Length for finding firing rate
NTrain=round(length(tt)*1); % number of training samples
NTest=round(length(tt)*0);% number of test samples
Trial_num=TotCells;

%% Test and Train data 
TrainInt=1:NTrain/LenWind;
% TesInt=NTrain/LenWind+1:NTrain/LenWind+NTest/LenWind;
TesInt=TrainInt;
yEstimated=zeros(length(TesInt),1);
TesIntTime=TesInt.*(tt(2)-tt(1));
%% Calculating firing rate and preprocessing the data
[ySum,yMean,yCells]=MeanFiringRate(F_binary, LenWind); % find mean of firing rate total
yMean=(yMean-min(yMean))/(max(yMean)-min(yMean));

[ySumSync,yMeanSync,yCellsSync]=MeanFiringRate(F_binarySync, LenWind); % find mean of firing rate Sync
yMeanSync=(yMeanSync-min(yMean))/(max(yMean)-min(yMean));

[ySumAsync,yMeanAsync,yCellsAsync]=MeanFiringRate(F_binaryAsync, LenWind); % find mean of firing rate Asunc
 yMeanAsync=(yMeanAsync-min(yMean))/(max(yMean)-min(yMean));
 
 Signal_M=conv(Signal_Mix,ones(LenWind,1)/LenWind,'same'); % same length 
 Signal_M=downsample(Signal_M,LenWind);
 minMix=min(Signal_M);
 maxMix=max(Signal_M);
 Signal_M=(Signal_M-minMix)/(maxMix-minMix);
 

 Signal_S=conv(Signal_slow(1:end-1)+midpoint,ones(LenWind,1)/LenWind,'same');
  minS=min(Signal_S);
 maxS=max(Signal_S);
 Signal_S=downsample(Signal_S,LenWind);
 Signal_S=(Signal_S-minS)/(maxS-minS);
 

 Signal_F=conv(Signal_fast(1:end-1),ones(LenWind,1)/LenWind,'same');
  minF=min(Signal_F);
 maxF=max(Signal_F);
 Signal_F=downsample(Signal_F,LenWind);
 Signal_F=(Signal_F- minF)/(maxF- minF);


%% Fit Kalman Filter for Mixed Signal
X=Signal_M(TrainInt)'; % is K * 1
Z=yMean(TrainInt,:)'; % is K * TotCells
%number of time bins
nt=length(X);
%Calculate the transition matrix (from x_t to x_t+1) using least-squares, and compute its covariance
X2 = X(2:end);
X1 = X(1:end-1);
A=X2*X1'*inv(X1*X1'); %Transition matrix
W=(X2-A*X1)*(X2-A*X1)'/(nt-1) ;%Covariance of transition matrix. Note we divide by nt-1 since only nt-1 points were used in the computation (that's the length of X1 and X2). We also introduce the extra parameter C here.
% Calculate the measurement matrix (from x_t to z_t) using least-squares, and compute its covariance
% In our case, this is the transformation from kinematics to spikes
H = Z*X'*(inv(X*X')); %Measurement matrix
Q = ((Z - H*X)*((Z - H*X)')) / nt; %Covariance of measurement matrix
%% Decode Mix spike by the KF
X=Signal_M(TesInt)';
Z=yMean(TesInt,:)';
%Initializations
num_states=1; %Dimensionality of the state
states=zeros(size(X)); %Keep track of states over time (states is what will be returned as y_test_predicted)
P_m=zeros([num_states,num_states]);
P=zeros([num_states,num_states]);
state=X(1); %Initial state
%Get predicted state for every time bin
        for t= 1:length(X)-1
            %Do first part of state update - based on transition matrix
            P_m=A*P*A'+W;
            state_m=A*state;

            %Do second part of state update - based on measurement matrix
            K=P_m*H'*inv(H*P_m*H'+Q); %Calculate Kalman gain
            P=(eye(num_states)-K*H)*P_m;
            state=state_m+K*(Z(t+1)-H*state_m);
            states(t+1)=squeeze(state); %Record state at the timestep
        end
y_test_predicted_Mixed=states';


%% Decode Sync spike by the KF
X=Signal_M(TesInt)';
Z=yMeanSync(TesInt,:)';
%Initializations
num_states=1; %Dimensionality of the state
states=zeros(size(X)); %Keep track of states over time (states is what will be returned as y_test_predicted)
P_m=zeros([num_states,num_states]);
P=zeros([num_states,num_states]);
state=X(1); %Initial state
%Get predicted state for every time bin
        for t= 1:length(X)-1
            %Do first part of state update - based on transition matrix
            P_m=A*P*A'+W;
            state_m=A*state;

            %Do second part of state update - based on measurement matrix
            K=P_m*H'*inv(H*P_m*H'+Q); %Calculate Kalman gain
            P=(eye(num_states)-K*H)*P_m;
            state=state_m+K*(Z(t+1)-H*state_m);
            states(t+1)=squeeze(state); %Record state at the timestep
        end
y_test_predicted_Sync=states';

%% Decode ASync spike by the KF
X=Signal_M(TesInt)';
Z=yMeanAsync(TesInt,:)';
%Initializations
num_states=1; %Dimensionality of the state
states=zeros(size(X)); %Keep track of states over time (states is what will be returned as y_test_predicted)
P_m=zeros([num_states,num_states]);
P=zeros([num_states,num_states]);
state=X(1); %Initial state
%Get predicted state for every time bin
        for t= 1:length(X)-1
            %Do first part of state update - based on transition matrix
            P_m=A*P*A'+W;
            state_m=A*state;

            %Do second part of state update - based on measurement matrix
            K=P_m*H'*inv(H*P_m*H'+Q); %Calculate Kalman gain
            P=(eye(num_states)-K*H)*P_m;
            state=state_m+K*(Z(t+1)-H*state_m);
            states(t+1)=squeeze(state); %Record state at the timestep
        end
y_test_predicted_Async=states';
%% Denormalize
% fprintf('Performance of the decoder model, mixedSignal, R2 =%f \n',rsquare(y_test_predicted_Mixed,Signal_M));
% fprintf('Performance of the decoder model, SlowSignal, R2 =%f \n',rsquare(y_test_predicted_Async,Signal_S));
% fprintf('Performance of the decoder model, FastSignal, R2 =%f \n',rsquare(y_test_predicted_Sync,Signal_F));
y_test_predicted_Mixed=(y_test_predicted_Mixed)*(maxMix-minMix)+minMix;
Signal_M=(Signal_M)*(maxMix-minMix)+minMix;

Signal_S=(Signal_S)*(maxF-minS)+minS;
y_test_predicted_Async=(y_test_predicted_Async)*(maxS-minS)+minS;

Signal_F=(Signal_F)*(maxF-minF)+minF;
y_test_predicted_Sync=(y_test_predicted_Sync)*(maxF-minF)+minF;

%% plot the result
figure
title('decode from Mixed Spikes','fontsize',22)
hold on
plot(tt(TesInt),Signal_M(TesInt),'linewidth',2);
plot(tt(TesInt),y_test_predicted_Mixed, ':r','linewidth',2);
hold off
legend({'true','estimated'},'fontsize',16);
xlabel('time (ms)','fontsize',22);
ylim([minMix,maxMix]);
set (gca, 'fontsize', 15)
figure
title('decode from Sync Spikes','fontsize',22)
hold on
plot(tt(TesInt),Signal_F(TesInt),'linewidth',2);
plot(tt(TesInt),y_test_predicted_Sync, '-.r','linewidth',2);
hold off
xlabel('time (ms)','fontsize',22);
ylim([minF,maxF]);
set (gca, 'fontsize', 15)
figure
title('decode from Async Spikes','fontsize',22)
hold on
plot(tt(TesInt),Signal_S(TesInt),'linewidth',2);
plot(tt(TesInt),y_test_predicted_Async, '-.r','linewidth',2);
hold off
xlabel('time (ms)','fontsize',22);
ylim([minS,maxS]);
set (gca, 'fontsize', 15)
% %% Distribution Estimation of spikes
% figure
% %%%%%%%%%%%%% Sync %%%%%%%%%
% subplot(511)
% syncDist='poisson'; % distribution type
% xSync=min(ySumSync):.1:max(ySumSync);
% hSync=hist(ySumSync,length(xSync)); % actual distribution
% hSync=hSync./max(hSync);
% dSync=fitdist(ySumSync,syncDist ) % parametric estimation of the distribution
% ySync = pdf(dSync,xSync);
% hold on % compare the results
% plot(xSync,hSync,'-b','LineWidth',5)
% plot(xSync,ySync,'-r','LineWidth',2)
% legend('true','estimated');
% hold off
% R2=squre(ySync(1:end),hSync(1:end));
% title(['Distribution Sync Spikes, dist=',syncDist,' with lambda=', num2str(dSync.lambda), ' , R2=', num2str(R2)])
% %%%%%%% Sync2 %%%%%%
% 
% subplot(512)
% syncDist='binomial';
% xSync=min(ySumSync):.1:max(ySumSync);
% hSync=hist(ySumSync,length(xSync));
% hSync=hSync./max(hSync);
% [phatSync, pciSync]=binofit(ySumSync,max(ySumSync)*ones(size(ySumSync)) );
% ySync = binopdf(xSync,mean(max(ySumSync)*ones(size(xSync))),mean(phatSync));
% hold on
% plot(xSync,hSync,'-b','LineWidth',5)
% plot(xSync,ySync,'-r','LineWidth',2)
% legend('true','estimated');
% hold off
% R2=squre(ySync(1:end),hSync(1:end));
% title(['Distribution Sync Spikes, dist=',syncDist,' with p=', num2str(mean(phatSync)), ' , R2=', num2str(R2)])
% 
% %%%%%%%%%%%%% Async %%%%%%%%%
% subplot(513)
% asyncDist='poisson';
% xAsync=min(ySumAsync):.1:max(ySumAsync);
% hAsync=hist(ySumAsync,length(xAsync));
% hAsync=hAsync./max(hAsync);
% dAsync=fitdist(ySumAsync,asyncDist )
% yAsync = pdf(dAsync,xAsync);
% hold on
% plot(xAsync,hAsync,'-b','LineWidth',5)
% plot(xAsync,yAsync,'-r','LineWidth',2)
% legend('true','estimated');
% hold off
% R2=squre(yAsync,hAsync);
% title(['Distribution Async Spikes, dist=',asyncDist,' with lambda=', num2str(dAsync.lambda), ' , R2=', num2str(R2)])
% % estimate async
% subplot(514)
% asyncDist='normal';
% xAsync=min(ySumAsync):.1:max(ySumAsync);
% hAsync=hist(ySumAsync,length(xAsync));
% hAsync=hAsync./max(hAsync);
% dAsync=fitdist(ySumAsync,asyncDist )
% yAsync = pdf(dAsync,xAsync);
% hold on
% plot(xAsync,hAsync,'-b','LineWidth',5)
% plot(xAsync,yAsync,'-r','LineWidth',2)
% legend('true','estimated');
% hold off
% R2=squre(yAsync,hAsync);
% title(['Distribution Async Spikes, dist=',asyncDist,' with sigma=', num2str(dAsync.sigma),' and mu=',num2str(dAsync.sigma) ' , R2=', num2str(R2)])
% 
% 
% %%%%%%%%%%%%% Mixed %%%%%%%%%
% subplot(515)
% aMixDist='poisson';
% x=min(ySum):.1:max(ySum);
% hMix=hist(ySum,length(x));
% hMix=hMix./max(hMix);
% dMix=fitdist(ySum,aMixDist )
% y = pdf(dMix,x);
% hold on
% plot(x,hMix,'-b','LineWidth',5)
% plot(x,y,'-r','LineWidth',2)
% legend('true','estimated');
% hold off
% R2=squre(y,hMix);
% title(['Distribution Mixed Spikes, dist=',aMixDist,' with lambda=', num2str(dMix.lambda), ' , R2=', num2str(R2)])
% 

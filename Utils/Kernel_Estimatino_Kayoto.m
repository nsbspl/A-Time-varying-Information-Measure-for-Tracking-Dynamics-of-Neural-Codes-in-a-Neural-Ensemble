%% TUTORIAL OF KERNEL DENSITY ESTIMATION BY SSKERNEL.M AND SSVKERNEL.M
clear all; close all;clc
% load('Run5_Data.mat')
%% Create samples
% First, let's create a sample data set. 
% Here `x' is a vector containing N samples, following a specified
% distribution. 
dt=tt(2)-tt(1);
LenWind=500;
TotCells=1;
%% Calculating firing rate and preprocessing the data
[ySum,yMean,yCells]=MeanFiringRate(F_binary, LenWind); % find mean of firing rate total
yMean=resample(yMean, LenWind*10,10);
x=tt(find(sum(F_binary,2) ~=0));
yMean=(yMean)/(max(yMean));


%% Plot a histogram
% Let's plot the data using a histogram. Here to see the noisy data 
% data structure, we use the bin-width 5 times smaller than an optimal 
% bin-width selected by histogram optimization method, 'sshist(x)'. 
subplot(3,1,1:2); hold on; 
plot(tt,yMean,'linewidth',2,'color',[1,0,1]);
ylim([0,1]);

%% Create a vector of estimation points
% Now, we estimate the underlying density from the samples.  
t = tt(1:LenWind:end);        %points at which density estimation is made
                            %points must be equi-distant.

%% Kernel Density Estimation by a Fixed Bandwidth
% To obtain a kernel estimated density with a fixed kernel bandwidth, type
[yf,tf,optw,w,c] = sskernel(x,t);

% yf is the estimated density with the optimal bandwidth optw. The
% estimation was given at the points in tf, which was automatically
% determined from the data. 
yf=(yf-min(yf))/(max(yf)-min(yf));
% Display the results.
subplot(3,1,1:2);
hold on; plot(tf,yf,':b','LineWidth',2);
set(gca,'XLim',[min(t) max(t)]);
ylabel('density');
subplot(3,1,3); hold on; 
plot(t,optw*ones(1,length(t)),':b','LineWidth',2);
set(gca,'XLim',[min(t) max(t)]);
ylabel('bandwidth');
drawnow;

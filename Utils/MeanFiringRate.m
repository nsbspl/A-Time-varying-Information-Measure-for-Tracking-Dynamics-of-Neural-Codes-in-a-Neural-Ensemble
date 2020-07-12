function [ySum,yMean, yCells ]= R_MeanFiringRate(Spike_Binary, TimeWindow)
% calculate the mean firing rate of spikes
% Spike_Binary is spiking activity
% TimeWindow smoothing window


% yCells firing rate each neuron
% yMean mean of firing rate
% ySum sum of all activities

yCells=conv2(Spike_Binary,ones(TimeWindow,1),'same');
yCells=downsample(yCells,TimeWindow);
yMean=mean(yCells,2);
ySum=sum(yCells,2);
%X=conv(X,ones(LenWind,1)/LenWind,'same');
%X=downsample(X,LenWind);
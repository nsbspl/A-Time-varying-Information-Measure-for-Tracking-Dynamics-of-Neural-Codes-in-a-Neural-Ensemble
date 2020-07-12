function [ySum,yMean, yCells ]= smoothNeuralActivity_v2(Spike_Binary, TimeWindow)
[nt,nc]=size(Spike_Binary);
for i=1:nc
yCells(:,i)=conv(Spike_Binary(:,i),TimeWindow,'same');

end
yMean=mean(yCells,2);
ySum=sum(Spike_Binary,2);
%X=conv(X,ones(LenWind,1)/LenWind,'same');
%X=downsample(X,LenWind);
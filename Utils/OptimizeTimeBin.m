function [optN,C,N]=R_OptimizeTimeBin(x,N,n,dt)
% Minimum number of bins (integer)
% N_MIN must be more than 1 (N_MIN > 1).
% Maximum number of bins (integer)
x_min = min(x);
x_max = max(x);

                   % # of Bins


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the Cost Function
for i = 1: length(N)
  
	xb=reshape(x(1:floor(length(x)/N(i))*N(i)),N(i),floor(length(x)/N(i)));
	ki = sum(xb,2);            % Count # of events in bins
	
	k = mean(ki);                   % Mean of event count
	v = sum( (ki-k).^2 )/N(i);      % Variance of event count

	C(i) = ( 2*k - v ) / (n*dt*(N(i))).^2;    % The Cost Function

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimal Bin Size Selectioin
[Cmin idx] = min(C);
optN = N(idx);  
% edges = linspace(x_min,x_max,N(idx)+1);  % Optimal segmentation
% % Display an Optimal Histogram and the Cost Function
% figure
% subplot(1,2,1); hist(x,edges); axis square;
% subplot(1,2,2); plot(D,C,'k.',optD,Cmin,'r*'); axis square;
% 

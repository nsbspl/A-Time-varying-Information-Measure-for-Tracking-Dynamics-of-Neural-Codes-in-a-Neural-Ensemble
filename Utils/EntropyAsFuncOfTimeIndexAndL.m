close all 
clear all
clc
%% Create or Load Simulated Stimulus
% load('Run4_Data.mat')
%%

WordLengths=[ 8,5 ]; % word lengths
pdt=(tt(2)-tt(1)); % calculate dt
%% Calculate Entropy as function of  
% window with length L 
% (dt is constant here and is the maximium resolution we can use)
Entropy={}
p=0;
for i =WordLengths
    p=p+1
    Entropy{p} = entropyFuncTimeIndexWindoL(F_binary,i,pdt);
end
%% Pre Plot (make the graphs smoother)
p=0;
for i =WordLengths
    p=p+1;
    Entropy{p} =resample(Entropy{p}, i*10,10);
end
%% Plot the result
subplot(length(WordLengths)+1,1,1);
 plot(tt.*pdt,Signal_Mix);
 ylim([min(Signal_Mix),max(Signal_Mix)]);
 title('signal Mixed');
 for i = 1 : length(WordLengths)
  subplot(length(WordLengths)+1,1,i+1);
 plot(tt.*pdt,Entropy{i});
 ylim([0,max(Entropy{i})]);
 title(['Entropy for word Length', num2str(WordLengths(i))]);
 end
 xlabel('time (ms)');
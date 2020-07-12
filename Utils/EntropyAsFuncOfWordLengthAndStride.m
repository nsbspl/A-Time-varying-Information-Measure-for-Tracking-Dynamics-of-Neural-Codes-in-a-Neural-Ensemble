close all 
clear all
clc
%% Create or Load Simulated Stimulus
% load('Run4_Data.mat')
%%
WordLengths=[ 5 ]; % word lengths
pdt=(tt(2)-tt(1)); % calculate dt
%% Calculate Entropy as function of  
% window with length L 
% (dt is constant here and is the maximium resolution we can use)
Entropy={}
p=0;
for i =WordLengths
    p=p+1
    Entropy{p} = entropyFuncStrideAndWindoL(F_binary,i,pdt);
end

%% Plot the result
subplot((WordLengths)+1,1,1);
 plot(tt.*pdt,Signal_Mix);
 ylim([min(Signal_Mix),max(Signal_Mix)]);
 title('signal Mixed');
 Res=Entropy{1};
 for i = 1 : (WordLengths)
 subplot((WordLengths)+1,1,i+1);
 plot(tt.*pdt,Res(:,i));
%  ylim([0,max(Res(:,i))]);
 title(['Entropy for Stride', num2str(i)]);
 end
 xlabel('time (ms)');
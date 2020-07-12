close all 
clear all
clc
%% Create or Load Simulated Stimulus
% load('Run4_Data.mat')
%%
WordLengths=[ 10, 5, 2 ]; % word lengths
pdt=(tt(2)-tt(1)); % calculate dt
DeltaT=[.05, .2, .5, .8, 1, 4, 8, 20, 50, 100, 200, 400, 500]; % temporal resolutions
%% Temporal resolution
[signalLength, NumCells] = size(F_binary); % signalLength number of samples, NumCells number of neurons
q=0;
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
Entropy={}
p=0;
for i =WordLengths
    p=p+1
    Entropy{p} = entropyFuncTimeIndexWindoL(Sp,i,dt);
end
%% Pre Plot (make the graphs smoother)
p=0;
for i =WordLengths
    p=p+1;
    Entropy{p} =resample(Entropy{p}, BinLength*i*10,10);
end
%% Plot the result
figure

subplot(length(WordLengths)+1,1,1);
 plot(Signal_Mix);
 ylim([min(Signal_Mix),max(Signal_Mix)]);
 title('signal Mixed');
 for i = 1 : length(WordLengths)
  subplot(length(WordLengths)+1,1,i+1);
 plot(Entropy{i});
%  ylim([0,max(Entropy{i})]);
 title(['Entropy for word Length=', num2str(WordLengths(i)),',timeRes=',num2str(dt)]);
 end
 xlabel('time (ms)');
end

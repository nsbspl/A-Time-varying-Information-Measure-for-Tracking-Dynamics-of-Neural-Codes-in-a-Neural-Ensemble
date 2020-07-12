close all 
clear all
clc
%% Create or Load Simulated Stimulus
load('Run4_Data.mat')
%%
WordLengths=[ 20,10, 5, 2,1 ]; % word lengths
pdt=(tt(2)-tt(1)); % calculate dt
DeltaT=[.05, .1, .2, .5, .8, 1, 4, 8, 20, 50, 100, 200]; % temporal resolutions
startInd=1;
endInd=length(F_binary);
%% Temporal resolution
Signal_fast=Signal_fast(1:end-1,:);
Signal_slow=Signal_slow(1:end-1,:);
[signalLength, NumCells] = size(F_binary(startInd:endInd,:));
normSignal_Mix=(Signal_Mix(startInd:endInd,:)-min(Signal_Mix(startInd:endInd,:)))/(max(Signal_Mix(startInd:endInd,:))-min(Signal_Mix(startInd:endInd,:)));
normSignal_fast=(Signal_fast(startInd:endInd,:)-min(Signal_fast(startInd:endInd,:)))/(max(Signal_fast(startInd:endInd,:))-min(Signal_fast(startInd:endInd,:)));
normSignal_slow=(Signal_slow(startInd:endInd,:)-min(Signal_slow(startInd:endInd,:)))/(max(Signal_slow(startInd:endInd,:))-min(Signal_slow(startInd:endInd,:)));
corResult=zeros(length(DeltaT),length(WordLengths),2);
q=0;
Entropy={}
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
    p=p+1
    Entropy{q,p} = entropyFuncTimeIndexWindoL(Sp,i,dt);
end

end
%%
for p=1:length(WordLengths)
    for q=1:length(DeltaT)
        BinLength = (DeltaT(q)/pdt);
         Entropy2{q,p} =resample(Entropy{q,p}, BinLength*WordLengths(p)*10,10);
    end
end
%% calculate Total entropy
load('Run4_Data_noise.mat')
F_binary_noise=F_binary_noise(1:end-1,:);
WordLengths_t=[10,1 ]; % word lengths
pdt=(tt(2)-tt(1)); % calculate dt
DeltaT_t=[.05, 20]; % temporal resolutions
q=0;
Entropy_t={}
for dt = DeltaT_t
    q=q+1;
    BinLength = (dt/pdt);
    Sp=zeros(signalLength / BinLength,NumCells);
    for i = 1: signalLength / BinLength
        for p=1:NumCells
            if sum( F_binary_noise( (i-1)*BinLength+1: (i)*BinLength,p)) >=1
                Sp(i,p) =1;
            end
        end
    end
    %% Calculate Entropy as function of  
% window with length L 
% (dt is constant here and is the maximium resolution we can use)

p=0;
for i =WordLengths_t
    p=p+1
    Entropy_t{q,p} = entropyFuncTimeIndexWindoL(Sp,i,dt);
end

end
%%

%% calculate sync info
infoSync=Entropy_t{1,1}-Entropy{1,2};
fprintf('maximum information of Sync = %d \n',mean(infoSync));
% calculate Async info
infoAsync=Entropy_t{2,2}-Entropy{9,5};
fprintf('maximum information of Async = %d \n',mean(infoAsync));
%%
En=[];
for p=1:length(WordLengths)
    for q=1:length(DeltaT)
     En(q,p)=mean(Entropy{q,p});
   
    end
end
figure

imagesc(1./WordLengths, log10(DeltaT), En)
xlabel('1/L');
ylabel( 'log(dt)');
title('Mean of Entropy');
ax = gca; % current axes
ax.FontSize = 15;
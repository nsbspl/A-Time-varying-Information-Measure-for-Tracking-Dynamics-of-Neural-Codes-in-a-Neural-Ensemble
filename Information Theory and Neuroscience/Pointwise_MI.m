function [MI, Ent] = Pointwise_MI (True_Event, Est_Event,p)
TW_bin = 5; %  msec resolution
bin_size = TW_bin/p.dt;
bin_range=1:bin_size:length(True_Event);
%% Find all spikes withing the range
T_event = TW_bin;
kernel_rectangle = [zeros(1/p.dt,1); ones(T_event/p.dt,1); zeros(1/p.dt,1)]; L_k = length(kernel_rectangle);

sig_True = zeros(size(True_Event));
sigg_ = conv(True_Event,kernel_rectangle); sig_True = sigg_(L_k/2:end-L_k/2);

sigg_ = [];
sig_Est = zeros(size(Est_Event));
sigg_ = conv(Est_Event,kernel_rectangle); sig_Est = sigg_(L_k/2:end-L_k/2);

count_true = sum(True_Event);
count_total = sum(Est_Event);
gg = diff(sig_Est.*sig_True);
TP = length(find(gg==1));
FN = count_true - TP; % total number of true events - true estimated events
FP = count_total - TP;
%%

Q=FP; % False Positive
R=FN; % False negative
S=TP; % True positive
P=length(bin_range) - (FP + FN + TP); % True negatives
       if P==0;
           P1=1; 
        else
           P1=P;
        end
       
        if Q==0;
            Q1=1; 
        else
            Q1=Q;
        end
        
        if R==0;
            R1=1; 
        else
            R1=R;
        end
        
        if S==0;
            S1=1; 
        else
            S1=S;
        end
        
 J = [P1,Q1;R1,S1]/length(bin_range);
 MI = sum(sum(J.*log2(J./(sum(J,2)*sum(J,1))))) * 1e3/TW_bin;
 
p1 = count_true/length(bin_range);
p0 = 1-p1;
Ent = - (p1*log2(p1)+p0*log2(p0)) * 1e3/TW_bin;


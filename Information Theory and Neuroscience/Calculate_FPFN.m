function [count_total,FP,FN,count_true,TP] = Calculate_FPFN(out_sync,INP,window_event,p)
%indx_,T_sp
T_event = window_event;
True_event = zeros(size(out_sync));
[i,j] = findpeaks(INP);
True_event(j) = 1;

kernel_rectangle = [zeros(1/p.dt,1); ones(T_event/p.dt,1); zeros(1/p.dt,1)]; L_k = length(kernel_rectangle);

sig_True = zeros(size(out_sync));
sigg_ = conv(True_event,kernel_rectangle); sig_True = sigg_(L_k/2:end-L_k/2);

sigg_ = [];
sig_Est = zeros(size(out_sync));
sigg_ = conv(out_sync,kernel_rectangle); sig_Est = sigg_(L_k/2:end-L_k/2);

count_true = sum(True_event);
count_total = sum(out_sync);
gg = diff(sig_Est.*sig_True);
TP = length(find(gg==1));
FN = count_true - TP; % total number of true events - true estimated events
FP = count_total - TP; % total number of estimated events - true estimated events
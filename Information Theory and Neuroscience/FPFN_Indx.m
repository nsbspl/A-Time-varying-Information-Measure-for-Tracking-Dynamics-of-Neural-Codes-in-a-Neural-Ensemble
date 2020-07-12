function [indx_TP,indx_FP,indx_FN] = FPFN_Indx(out_sync,INP,window_event,p)
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


gg_TP = (out_sync.*sig_True);%diff(sig_Est.*sig_True);%
indx_TP = find(gg_TP>=1);% + T_event/p.dt/2;%find(gg_TP==1);% %to find the center of the event
gg_FP = out_sync.*(~sig_True);
indx_FP = find(gg_FP>=1);
gg_FN = (~sig_Est).*INP;
indx_FN = find(gg_FN>=1);
%%% Remove events whose ISIs < 50 msec
if (length(indx_TP) == 0)
   indx_TP = [];
else
    indx_no = [indx_TP(1);find(diff(indx_TP) < 50/p.dt)];
    indx_TT = indx_TP; indx_TT(indx_no) = 0;
    indx_TP = []; indx_TP = indx_TT(indx_TT>0);
end




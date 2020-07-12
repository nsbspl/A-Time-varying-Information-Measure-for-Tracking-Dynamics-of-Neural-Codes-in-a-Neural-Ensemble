function [F_mix, F_slow, F_fast] = Reconstruct_2_window(psth_tot, Filters_, out_sync, window_event,p)
q =2;
L_W = floor(length(Filters_)/2);
F_ = zeros(length(Filters_)+length(psth_tot)-1,q);
T_event = window_event;
% Rectangular Window
kernel_rectangle = [zeros(1/p.dt,1); ones(T_event/p.dt,1); zeros(1/p.dt,1)]; L_k = length(kernel_rectangle);
sigg_ = [];
sig_Est = zeros(size(psth_tot));
sigg_ = conv(out_sync,kernel_rectangle); sig_Est = sigg_(L_k/2:end-L_k/2);
Sig_indx = zeros(length(psth_tot),q);
Sig_indx(:,1) = ~sig_Est;
Sig_indx(:,2) = sig_Est;
for i=1:q
    aa_ = zeros(size(psth_tot));
    %ind_flag = isempty(ind{i});
    aa_= psth_tot.*Sig_indx(:,i);
    F_(:,i) = conv(aa_-mean(aa_),Filters_(:,i));
end
FF = sum(F_,2);
F_mix = FF(L_W+1:end-L_W);
F_mix = F_mix(1:length(psth_tot));
F_slow = F_(L_W+1:end-L_W,1);
F_slow = F_slow(1:length(psth_tot));
F_fast = F_(L_W+1:end-L_W,2);
F_fast = F_fast(1:length(psth_tot));
    
    
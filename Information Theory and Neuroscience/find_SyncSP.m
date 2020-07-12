function [indx_async,indx_sync,sync_event,M_S] = find_SyncSP(psth_Learning,percentage_sync,window_sync,TW,Trial_num,p)
%indx_,T_sp
minNumSP = percentage_sync*Trial_num;
T_sync = window_sync;
%% rectangular PSTH
kernel_rectangle = [zeros(1/p.dt,1); ones(T_sync/p.dt,1); zeros(1/p.dt,1)]; L_k = length(kernel_rectangle);
% first we make a simple PSTH with window size of window_sync msec.
% Note that this PSTH needs to be shifted appropriately. So for this
% reason, we make the Kernel PSTH as well.
sig_PSTH = zeros(size(psth_Learning));
sigg_ = conv(psth_Learning,kernel_rectangle); sig_PSTH = sigg_(L_k/2:end-L_k/2);
%% Gaussian Kernel PSTH
FG_K_PSTH = zeros(size(psth_Learning));
Tw_new = TW;% in msec   2*1e1*p.dt;% Tw is in ms
bin_width = Tw_new;%/p
res = p.dt;
t=-5*bin_width:res:5*bin_width;
Gauss_kernel = 1/sqrt(2*pi*(bin_width)^2)*exp(-t.^2/(2*(bin_width)^2)); 
Gauss_kernel = Gauss_kernel/(norm(Gauss_kernel,1)*res);
L_Hkernel = length(Gauss_kernel);
Kernel_psth_new = conv(psth_Learning,Gauss_kernel);%/(N*res*1e-3); % KHz
Kernel_psth_new = Kernel_psth_new(L_Hkernel/2:end-L_Hkernel/2);
FG_K_PSTH = Kernel_psth_new(1:length(sig_PSTH));
%% Correct the shift
FG_PSTH = zeros(size(psth_Learning));
xc = xcorr(sig_PSTH/std(sig_PSTH),FG_K_PSTH/std(FG_K_PSTH),200/p.dt);
cnt = 200/p.dt;
[val,ind] = max(xc); 
L2 = length(sig_PSTH) - abs(ind-cnt) + 1;
FG_PSTH(1:L2) = sig_PSTH(abs(ind-cnt):end);
%% 
FG_New = FG_K_PSTH/max(FG_K_PSTH) * max(FG_PSTH); % FG_PSTH
[i,j] = findpeaks(FG_New);
%%% --- Note 
sync_event = j(i>=minNumSP);
% sync_event = j(i>=percentage_sync*max(FG_PSTH));

MM_S = zeros(size(FG_PSTH));
MM_S(sync_event) = 1;
pf = conv(MM_S,kernel_rectangle); M_S = pf(L_k/2:end-L_k/2); % This includes all the time around the peak of the sync_event
indx_sync = find(psth_Learning.*M_S == 1);
indx_async = find((psth_Learning - psth_Learning.*M_S) == 1);

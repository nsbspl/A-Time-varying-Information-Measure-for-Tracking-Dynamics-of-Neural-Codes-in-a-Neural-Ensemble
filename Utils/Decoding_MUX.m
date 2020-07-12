Signal_Inp = Signal_Mix;
%%%--- parameters pre-settings to separate synchronous and asynchronous spikes
Trial_num = 30; % Network size
percentage_sync = 0.40; % 40% of network size
window_sync = 5;% msec
window_event = 5; % msec
tw = 100; % STA time window(msec)
%% Calculating STA filter, (1) Sync- and Async- STAs & (2) All-spikes- STA
LearningTime = 7/p.dt *1e3; % 7 sec for training
INP_Learning = Signal_Inp(1:LearningTime);%
Inp_Test = Signal_Inp(LearningTime+1:end); s_orig = (Inp_Test - mean(Inp_Test)); %s_orig = s_orig/std(s_orig);
% read the spikes

F_binary;
num_trial = size(F_binary,2);        
%% Learning
TW = 1; % kernel-width (msec) for making instantaneous firing rate
psth_total = PSTH_(F_binary(:,1:num_trial), p.dt, p.dt); % superimpose spikes
[F_T, indx_T] = spike_times(F_binary(1:LearningTime,1:num_trial),p); % find spikes times to calculate STAs
psth_Learning = psth_total(1:LearningTime);
FG_mix = KernelPSTH (psth_Learning,TW,p.dt,num_trial);
TT = find(psth_Learning>0);
%--- find synchronous and asynchronous spikes
q=2;
dt = p.dt;
T_sp = [];
[indx_async,indx_sync,sync_event] = find_SyncSP(psth_Learning,percentage_sync,window_sync,TW,Trial_num,p);
T_sp{1} = indx_async;
T_sp{2} = indx_sync;
out_sync = zeros(length(psth_Learning),1);
out_sync(sync_event) = 1;  
[count_total,FP,FN,count_true] = Calculate_FPFN(out_sync,Fast_signal(1:LearningTime),window_event,p);
aa = find(psth_Learning>0);
L_tw = tw/p.dt;
%--- Calculate STAs for Synchronous and Asynchronous spikes
STA = zeros(2*tw/p.dt+1,q);
length_sta = 2*tw/p.dt+1;
for i_ = 1:q
spike_indx = T_sp{i_};
maxindx = aa(end);
indx_int = find(spike_indx>length_sta & spike_indx<maxindx-length_sta); a_ = isempty(indx_int);
if a_==1
    sta = zeros(length(STA),1);
    STA(:,i_) = sta;
else
    [sta, stc] = Spike_Triggered_Covariance2(T_sp{i_}*dt, aa(end)*dt, INP_Learning, tw, dt); %INP_Learning
    sta_ = sta - sta(1);
    STA(:,i_) = sta_/(norm(sta_,1)*p.dt);
end
end
% Find the best coefficients for each STA
sig_regul = INP_Learning - midpoint;
Inp_M = zeros(length(sig_regul),q);
for i_=1:q
    aa_ = zeros(size(psth_Learning));
    aa_(T_sp{i_}) = 1;
    Out_ = psth_Learning.*aa_;
    Inp_ = conv(STA(:,i_),(Out_ - mean(Out_)));%filter(sta,1,(Out_Test - mean(Out_Test)));
    Inp_ = (Inp_(1:length(Out_)))/Trial_num; %if norm(Inp_)~=0, Inp_ = Inp_Est_/std(Inp_Est_); end
    xc = xcorr(sig_regul/std(sig_regul),Inp_/std(Inp_),200/p.dt);
    cnt = 200/p.dt;
    [val,ind] = max(xc); 
    L2 = length(sig_regul) - abs(ind-cnt) + 1;
    Inp_M(1:L2,i_) = Inp_(abs(ind-cnt):end);
end
a_coef = Inp_M\sig_regul;
for i_=1:q
    STA(:,i_) = STA(:,i_)*a_coef(i_);
end
%--- Calculate STA for All spikes (conventional decoding)
[F_T, indx_T] = spike_times(F_binary(1:LearningTime,1:num_trial),p);
sp_Time = [];
for ik=1:num_trial
    k_ = [];
    k_ = find(F_T(:,ik)>0); 
    sp_Time = [sp_Time F_T(k_,ik)'];
end
sp_Time_L = [];
sp_Time_L = sort(sp_Time);
rate = p.dt;% least time window in msec
[sta_, stc] = Spike_Triggered_Covariance2(sp_Time_L, sp_Time_L(end), sig_regul, tw, p.dt);%/Trial_num;
sta = sta_ - sta_(1);
sta = sta/(norm(sta)*p.dt);
% modify the coefficient of the sta_total
Inp_sta = zeros(length(sig_regul),1);
Out_ = psth_Learning;
Inp_ = conv(sta,(Out_ - mean(Out_)));%filter(sta,1,(Out_Test - mean(Out_Test)));
Inp_ = (Inp_(1:length(Out_)))/Trial_num; %if norm(Inp_)~=0, Inp_ = Inp_Est_/std(Inp_Est_); end
xc = xcorr(sig_regul/std(sig_regul),Inp_/std(Inp_),200/p.dt);
cnt = 200/p.dt;
[val,ind] = max(xc); 
L2 = length(sig_regul) - abs(ind-cnt) + 1;
Inp_sta(1:L2) = Inp_(abs(ind-cnt):end);
a_sta = Inp_sta\sig_regul;
sta = sta*a_sta;

figure; hold on,
plot(-tw:p.dt:tw,STA)
legend('Asynch-STA','Synch-STA')
xlabel('Time (msec)')
ylabel('Amp (normalized)')
figure;
plot(-tw:p.dt:tw,sta,'k')
sta_conv = sta;
legend('All-Spike-STA')
xlabel('Time (msec)')
ylabel('Amp (normalized)')
%% Reconstruct Signals in Tested Simulation
psth_Test = psth_total(LearningTime+1:end);
T_sp = [];
[indx_async,indx_sync,sync_event] = find_SyncSP(psth_Test,percentage_sync,window_sync,TW,Trial_num,p);
T_sp{1} = indx_async;
T_sp{2} = indx_sync;
out_sync = zeros(length(psth_Test),1);
out_sync(sync_event) = 1;
%--- Reconstruct Input using Asynch & Synch STAs
[Inp_Est_M, Inp_slow, Inp_fast] = Reconstruct_2_window((psth_Test-mean(psth_Test)), STA, out_sync, window_event,p);
%--- Reconstruct Input using All-Spike-STA
L_tw = tw/p.dt;
Inp_sta = zeros(length(psth_Test),1);
Out_Test = psth_Test;
Inp_Est_ = conv(sta,(Out_Test - mean(Out_Test)));%filter(sta,1,(Out_Test - mean(Out_Test)));
Inp_sta = (Inp_Est_(L_tw+1:end-L_tw));% Inp_sta = Inp_sta/std(Inp_sta);
xc = xcorr(s_orig,Inp_sta,200/p.dt);
cnt = 200/p.dt;
[val,ind] = max(xc);
if abs(ind-cnt)~=0
L2 = length(s_orig) - abs(ind-cnt) + 1;
Inp_sta(1:L2) = Inp_sta(abs(ind-cnt):end); Inp_sta = Inp_sta/std(Inp_sta);
end
%% Compare the Reconstructed Inputs
figure; hold on,
tt_Test = 7e3:p.dt:p.tStop;
plot(tt_Test, s_orig/std(s_orig),'k')
plot(tt_Test,Inp_Est_M/std(Inp_Est_M))
plot(tt_Test,Inp_sta/std(Inp_sta),'r')
title('Original vs. Reconstructed Stimulus')
legend('Original','De-multiplexing','Conventioanl Decoding')
xlabel('Time (msec)')
ylabel('Amp (normalized)')
%% Raster plot of Spikes (with Synch- Asynch- spikes colored)
indx = (3e3/p.dt:5e3/p.dt); % 2 sec
figure; plot(p.dt:p.dt:length(indx)*p.dt,Signal_Mix(indx))
axis([0 length(indx)*p.dt -50 150]) 
xlabel('Time (msec)')
ylabel('Amp (pA)')
title('Stimulus')

psth_total = PSTH_(F_binary, p.dt, p.dt);
psth_ = psth_total(indx);
[indx_async,indx_sync,sync_event,M_S] = find_SyncSP(psth_total,percentage_sync,window_sync,TW,Trial_num,p);
sig_sync = M_S; 
sig_async = ~M_S;
Sig_T{1}=sig_async(indx);
Sig_T{2}=sig_sync(indx);
% 
% 
figure;
%ax(1) = subplot(2,1,1);
col = {'b','r'};%{'k','g','c','y','m','b'};
hold on,
for i=1:10%Trial_num
    for j=1:q
        jind = (find(F_binary(indx,i).*Sig_T{j}>0)*p.dt)';
        %plot(jind*p.dt,i*ones(length(jind),1),'.','Color',col{j})
        line(repmat(jind,2,1),repmat([i-1;i],1,length(jind)),'color',col{j})
    end
end
hold off
axis([0 length(indx)*p.dt 0 10]) 
xlabel('Time (msec)')
ylabel('Neuron Index')
title('Colored Raster-plot')
%% Firing rate
fr_neuron = sum(psth_total)/EndTime/Trial_num;
TW = 100; % wide time window (represents rate code) & TW = 5 (temporal code)
FG_mix = KernelPSTH (psth_total,TW,p.dt,Trial_num);
scale_ = fr_neuron/mean(FG_mix); % sum(scale_*FG_mix)*p.dt*1e-3/EndTime = sum(psth_total) / Trial_num
figure; plot(indx,scale_*FG_mix(indx),'k')
xlabel('Time (msec)')
ylabel('Firing Rate (Hz)')
title('Time-varying firing rate')

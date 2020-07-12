clear all
P1 = genpath('R:\MILAD\Results_March_2015\J_nCor_ModifiedSigma_40');% R:\MILAD\Results_March_2015\J_DynamicClamp
P2 = genpath('R:\MILAD\Results_March_2015\IEEE_CIM');%genpath('C:\Users\Milad Lankarany\Desktop\REcent\IEEE_CIM');
P3 = genpath('C:\Users\Milad Lankarany\Documents\MATLAB\Theoretical');
P4 = genpath('R:\MILAD\Results_March_2015\J_nCor_ModifiedSigma_40\2nd-March2016');
P5 = genpath('R:\MILAD\Results_March_2015\J_nCor_ModifiedSigma_40\3rd-March2016');
addpath(P1,P2,P3)
SamplingTime = 0.05;% msec
EndTime = 100;%sec
p = NeuronType(SamplingTime,EndTime,'CD');
tt = 0:p.dt:p.tStop;

midpoint = 15;
a_slow = 60;
a_fast = 85;
% midpoint = 70;
% a_slow = 120;%60;%
% a_fast = 170;%85;%

percentage_sync = 0.6;
window_sync = 5;% msec
%Trial_num = 10;
window_event = 5; % msec


WhichInp = 'mix';
    load Fast_signal, 
    load Slow_signal
    Signal_mix = (midpoint + a_slow*Slow_signal + a_fast*Fast_signal)*1e-3;
    Signal_Inp = Signal_mix;
    
%Sigma = (0:5:80)*3;%[0;15;60;120;180;240];%0:10:160;%[0;10;50;120];%[10;20;30;40;50;60;70;80;90;100;110];%[0.1; 0.2; 0.3; 0.4; 0.5];
Sigma = [0;0.25;0.5;0.75;1];
N = [1;2;5;10];%[1;2;5;10;20;50];
CC = [0.05; 0.1; 0.15; 0.2; 0.3; 0.4; 0.5];
tw = 100;
MI = zeros(2,length(Sigma),1);% 1st is for slow and second is for fast comparisons
FR = MI;
Err = MI;
MI_std = MI;
CF_mu_M = MI; CF_sem_M = MI;
CF_mu_S = MI; CF_sem_S = MI;
%%
LearningTime = EndTime/p.dt *1e3;
INP_Learning = [Slow_signal(1:LearningTime) Fast_signal(1:LearningTime)];%Signal_Inp(1:LearningTime);%
Inp_Test = Signal_Inp(LearningTime+1:end); s_orig = (Inp_Test - mean(Inp_Test)); s_orig = s_orig/std(s_orig);
Inp_Test_slow = Slow_signal(LearningTime+1:end); s_orig_slow = (Inp_Test_slow - mean(Inp_Test_slow)); s_orig_slow = s_orig_slow/std(s_orig_slow);
Inp_Test_fast = Fast_signal(LearningTime+1:end); s_orig_fast = (Inp_Test_fast - mean(Inp_Test_fast)); s_orig_fast = s_orig_fast/std(s_orig_fast);
%[f,Power_orig] = PSD(s_orig, p);
for i=length(N)%1:length(N)%
    Trial_num = N(i);
    for k = 4;%1:length(Sigma)%1:length(Sigma)%1:length(CC)
       %% Learning
        
        names = strcat('Fmix',{'_'},'sigma',num2str(k),{'_'},'N',num2str(i),{'.mat'}); % this is a cell and the 'save' only accpets char
        C_mix = char(names);
        load(C_mix); 
        F_binary = F_mix;
        
        names = strcat('Fmix',{'_'},'sigma',num2str(k),{'_'},'N',num2str(i),'_2nd',{'.mat'}); % this is a cell and the 'save' only accpets char
        C_mix = char(names);
        load(C_mix); 
        F_binary_2nd = PSTH_(F_mix, p.dt, p.dt);
        
        names = strcat('Fmix',{'_'},'sigma',num2str(k),{'_'},'N',num2str(i),'_3rd',{'.mat'}); % this is a cell and the 'save' only accpets char
        C_mix = char(names);
        load(C_mix); 
        F_binary_3rd = PSTH_(F_mix, p.dt, p.dt);
        
%%
TW = 1;
psth_total = PSTH_(F_binary, p.dt, p.dt);
[F_T, indx_T] = spike_times(F_binary(1:LearningTime,:),p);
psth_Learning = psth_total(1:LearningTime);
FG_mix = KernelPSTH (psth_Learning,TW,p.dt,Trial_num);
TT = find(psth_Learning>0);
% find different classes
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
%tw = 50; %msec
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
    [sta, stc] = Spike_Triggered_Covariance2(T_sp{i_}*dt, aa(end)*dt, INP_Learning(:,i_), tw, dt); %INP_Learning
    sta_ = sta - sta(1);
    STA(:,i_) = sta_/(norm(sta_,1)*p.dt);
end

end
% Find the best coefficients
sig_regul = Signal_mix(1:LearningTime) - mean(Signal_mix(1:LearningTime));%midpoint*1e-3;
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

% for All spikes
[F_T, indx_T] = spike_times(F_binary(1:LearningTime,:),p); % spike times are in msec
sp_Time = [];
for ik=1:Trial_num
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
%% Reconstruct Signals
if LearningTime == EndTime/p.dt *1e3;
    psth_Test = psth_total;
else
    psth_Test = psth_total(LearningTime+1:end);
end

T_sp = [];
[indx_async,indx_sync,sync_event] = find_SyncSP(psth_Test,percentage_sync,window_sync,TW,Trial_num,p);
T_sp{1} = indx_async;
T_sp{2} = indx_sync;
out_sync = zeros(length(psth_Test),1);
out_sync(sync_event) = 1;

Inp_Est_M = zeros(length(s_orig),q);
%Inp_Est_M = Reconstruct_2_modified((psth_Test-mean(psth_Test)),T_sp, STA); Inp_Est_M = Inp_Est_M/std(Inp_Est_M);
[Inp_Est_M, Inp_slow, Inp_fast] = Reconstruct_2_window((psth_Test-mean(psth_Test)), STA, out_sync, window_event,p);

L_tw = tw/p.dt;
Inp_sta = zeros(length(psth_Test),1);
Out_Test = psth_Test;
Inp_Est_ = conv(sta,(Out_Test - mean(Out_Test)));%filter(sta,1,(Out_Test - mean(Out_Test)));
Inp_sta = (Inp_Est_(L_tw+1:end-L_tw));% Inp_sta = Inp_sta/std(Inp_sta);

%% Calculate Mutual Information Lower Band
%Signal_mix = (midpoint + a_slow*Slow_signal + a_fast*Fast_signal)*1e-3;
[MI_sta,SNR_conv,f] = Mutual_Info_LB (zscore(Signal_mix), zscore(Inp_sta),p); 
MI_Mux = Mutual_Info_LB (zscore(Signal_mix), zscore(Inp_Est_M),p);
[MI_Async,SNR_Async,f] = Mutual_Info_LB (zscore(Slow_signal), zscore(Inp_slow),p);
MI_Sync = Mutual_Info_LB (Fast_signal/std(Fast_signal), Inp_fast/std(Inp_fast),p);
% [Inp_Est_M_slow, Inp_Est_M_fast] = separateSignal(Inp_Est_M, p);
% [Inp_sta_slow, Inp_sta_fast] = separateSignal(Inp_sta, p);
%% Calculate Mutual Information Upper Band
R_spk = zeros(length(psth_total),3);
R_spk(:,1) = psth_total; R_spk(:,2) = F_binary_2nd; R_spk(:,3) = F_binary_3rd;
[MI_UB, CRR_conv, f] = Mutual_Info_UB (R_spk,p);
% for Async spikes
R_spk_Async = R_spk;

T_sp = [];
[indx_async,indx_sync,sync_event] = find_SyncSP(R_spk(:,1),percentage_sync,window_sync,TW,Trial_num,p);
T_sp{1} = indx_async;
aa_ = zeros(size(psth_total));
aa_(T_sp{1}) = 1;
R_spk_Async(:,1) = R_spk(:,1).*aa_;

T_sp = [];
[indx_async,indx_sync,sync_event] = find_SyncSP(R_spk(:,2),percentage_sync,window_sync,TW,Trial_num,p);
T_sp{1} = indx_async;
aa_ = zeros(size(psth_total));
aa_(T_sp{1}) = 1;
R_spk_Async(:,2) = R_spk(:,2).*aa_;

T_sp = [];
[indx_async,indx_sync,sync_event] = find_SyncSP(R_spk(:,3),percentage_sync,window_sync,TW,Trial_num,p);
T_sp{1} = indx_async;
aa_ = zeros(size(psth_total));
aa_(T_sp{1}) = 1;
R_spk_Async(:,3) = R_spk(:,3).*aa_;

[MIasync_UB, CRR_Async, f] = Mutual_Info_UB (R_spk_Async,p);

%% CAlculate point-wise mutual information
T_event = window_event;
True_event = zeros(size(out_sync));
[i,j] = findpeaks(Fast_signal);
True_event(j) = 1;
[MI, H_entropy] = Pointwise_MI (True_event, out_sync,p);
    end
end

tx = 1:10*1e3/p.dt; % match with the data shown in Bayesian Decoding
figure;
plot(tx*p.dt,zscore(Signal_mix(tx)),'r','LineWidth',2)
hold on,
plot(tx*p.dt,zscore(Inp_sta(tx)),'k')

figure;
plot(tx*p.dt,zscore(Signal_mix(tx)),'r','LineWidth',2)
hold on,
plot(tx*p.dt,zscore(Inp_Est_M(tx)),'k')

%% /figures
Signal_mix = (midpoint + a_slow*Slow_signal + a_fast*Fast_signal)*1e-3;
cutoff1 = 0.1; cutoff2 = 1000; nfft = 1e6;

figure; % power spectrum slow signal
[Pslow,f] = pwelch(a_slow*Slow_signal,[],[],nfft,20e3); % pA
df = f(2) - f(1);
ni=(find(f<=cutoff2 & f>=cutoff1));
loglog(f(ni),df*Pslow(ni),'b') % power = pA^2/sec

figure; % power spectrum fast signal
[Pfast,f] = pwelch(a_fast*Fast_signal,[],[],nfft,20e3);
df = f(2) - f(1);
ni=(find(f<=cutoff2 & f>=cutoff1));
loglog(f(ni),df*Pfast(ni),'r')

figure; % power spectrum mixed signal
[Pmix,f] = pwelch(Signal_mix*1e3,[],[],nfft,20e3);
ni=(find(f<=cutoff2 & f>=cutoff1));
loglog(f(ni),df*Pmix(ni),'k')


figure; 
semilogx(f,log2(1+SNR_conv),'b')
hold on,
semilogx(f,-log2(1-sqrt(CRR_conv)),'r')

figure; 
semilogx(f,log2(1+SNR_Async),'b')
hold on,
semilogx(f,-log2(1-sqrt(CRR_Async)),'r')








clear all
%P1 = genpath('R:\MILAD\Results_March_2015\Info_Transmission\Slow_Signal (all-to-all conncetions- 5 sec)\Conductance'); %Slow_Signal (Sparse conncetions- 5 sec) Slow_Signal (Hybrid conncetions- 5 sec)    
%P1 = genpath('R:\MILAD\Results_March_2015\Info_Transmission\Mixed_Signal (all-to-all conncetions)\Conductance'); % Mixed signal
P1 = genpath('R:\MILAD\Results_March_2015\Info_Transmission\Mixed_Signal (all-to-all conncetions)\Saved (10L)- Slow'); % Slow
P2 = genpath('C:\Users\Milad Lankarany\Documents\MATLAB\Information Theory and Neuroscience');
P3 = genpath('R:\MILAD\Results_March_2015\IEEE_CIM');
addpath(P1,P2,P3)
SamplingTime = 0.05;% msec
EndTime = 5;%sec
p.dt = SamplingTime; %= NeuronType(SamplingTime,EndTime,'CD');
p.tStop = EndTime*1e3;
tt = 0:p.dt:p.tStop;
L = EndTime/p.dt*1e+3 + 1;
midpoint = 0;%15;
a_slow = 60;
a_fast = 0;% 85;
% midpoint = 70;  FOR EXPERIMENTAL DATA
% a_slow = 120;%60;%
% a_fast = 170;%85;%

percentage_sync = 0.30;
window_sync = 5;% msec
%Trial_num = 50;
window_event = 5; % msec
    load Slow_signal,
    Signal_mix = (midpoint + a_slow*Slow_signal(1:L))*1e-3;
    Signal_Inp = Signal_mix;

tau_rise = 0.5; tau_fall = 3; tEnd = 50*1e-3;
Gz = exp2syn_new(p.dt, tEnd, tau_rise, tau_fall); % tau_rise = 0.5 and tau_decay = 3;
L_Gz = length(Gz);  
    
Sigma = (0:5:40)*3;%(0:5:80)*3;%[0;15;60;120;180;240];%0:10:160;%[0;10;50;120];%[10;20;30;40;50;60;70;80;90;100;110];%[0.1; 0.2; 0.3; 0.4; 0.5];
Sigma(1) = 2*3;
N = [1;2;5;10;30;50;500];
tw = 100;
Num_layer = 4;
MI = zeros(Num_layer+1,6);%length(N)
FR = MI;
CF = MI;
FP = MI; TP = MI;
W_syn = zeros(N(end),Num_layer+1,6);
%%
LearningTime = 5/p.dt *1e3;
INP_Learning = Signal_Inp(1:LearningTime);%[Slow_signal(1:LearningTime) Fast_signal(1:LearningTime)];%
Inp_Test =  INP_Learning; s_orig = (Inp_Test - mean(Inp_Test)); %Signal_Inp(LearningTime+1:end); s_orig = s_orig/std(s_orig);

%%% I consider the PSTH of the layer 0  as the reference (instead of the
% stimulus) because no decoding algorithms is necessary to be used here.
i=length(N);
Trial_num = N(i);
load F_mix_layer_0
F_binary = F_mix;
TW_psth = 10;
psth_total = PSTH_(F_binary, p.dt, p.dt);
fr_neuron = sum(psth_total)/EndTime/Trial_num;
FG_mix = KernelPSTH (psth_total,TW_psth,p.dt,Trial_num);
scale_ = fr_neuron/mean(FG_mix);
FG_Ref = scale_*FG_mix;
%%%
s_orig = FG_Ref;

for mil = 1%0:Num_layer
    
for i=length(N)%1:length(N)%
    Trial_num = N(i);
    for k = 8%3:8%1:length(Sigma)%1:length(Sigma)%1:length(CC)
       %% Learning
        
       if mil==0
        load F_mix_layer_0
        F_binary = F_mix;
       else
        names = strcat('Fmix_2',{'_'},'sigma',num2str(k),{'_'},'layer',num2str(mil),{'.mat'});
        C_mix = char(names);
        load(C_mix); 
        F_binary = F_mix_2(:,1:Trial_num);
       end

%% Synaptic Weights
F_synaptic = zeros(size(F_binary));
      for kk=1:Trial_num     
          fr_fast = conv(F_binary(:,kk),Gz); fr_fast = fr_fast(L_Gz/2:end-L_Gz/2);
          F_synaptic(:,kk) = fr_fast;
      end
rand_spars = 1:Trial_num; % peak random number from [1-Trial_num]
a_scale = F_synaptic(:,rand_spars)\Signal_mix;
W_syn(:,mil+1,k) = a_scale;
%% Learning
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
% [count_total,fp,fn,count_true] = Calculate_FPFN(out_sync,Fast_signal(1:LearningTime),window_event,p);
FP (mil+1,k) = length(indx_sync);
% TP (mil,k) = count_true;

% aa = find(psth_Learning>0);
% % for All spikes
% [F_T, indx_T] = spike_times(F_binary(1:LearningTime,:),p); % spike times are in msec
% sp_Time = [];
% for ik=1:Trial_num
%     k_ = [];
%     k_ = find(F_T(:,ik)>0); 
%     sp_Time = [sp_Time F_T(k_,ik)'];
% end
% sp_Time_L = [];
% sp_Time_L = sort(sp_Time);
% rate = p.dt;% least time window in msec
% [sta_, stc] = Spike_Triggered_Covariance2(sp_Time_L, sp_Time_L(end), s_orig, tw, p.dt);%/Trial_num;
% sta = sta_ - sta_(1);
% sta = sta/(norm(sta)*p.dt);
% % modify the coefficient of the sta_total
% L_tw = tw/p.dt;
% Inp_sta = zeros(length(s_orig),1);
% Out_ = psth_Learning;
% Inp_ = conv(sta,(Out_ - mean(Out_)));%filter(sta,1,(Out_Test - mean(Out_Test)));
% Inp_sta = Inp_(L_tw+1:end-L_tw); %if norm(Inp_)~=0, Inp_ = Inp_Est_/std(Inp_Est_); end
% a_sta = Inp_sta\s_orig; %a_sta = std(sig_regul)/std(Inp_sta);
% sta = sta*a_sta; % Modified sta (conventional coding)
% Inp_sta = a_sta*Inp_sta;
%% Reconstruct Signals in Tested Simulation
fr_neuron = sum(psth_total)/EndTime/Trial_num;
FG_mix = KernelPSTH (psth_total,TW_psth,p.dt,Trial_num);
scale_ = fr_neuron/mean(FG_mix);
FG_test = scale_*FG_mix;
CF(mil+1,k) = (1 - norm(s_orig - FG_test)/norm(s_orig));%(1 - norm(s_orig - Inp_sta)/norm(s_orig));
FR(mil+1,k) = fr_neuron;
    end
end
end


layy = 1;
sig = 6;
a_dis = W_syn(:,layy+1,sig);
figure; hist(a_dis(a_dis>=0),20)
figure; hist(a_dis(a_dis<0),20)
%% For synaptic conductances
% ind_exc = find(a_scale>=0);
% ind_inh = find(a_scale<0);
% I_exc = F_synaptic(:,ind_exc)*a_scale(ind_exc);
% I_inh = F_synaptic(:,ind_inh)*a_scale(ind_inh);
% 
% figure; hold on,
% plot(inp_,'r')
% plot(I_exc + I_inh,'k--')
% 
% W_E = zeros(length(ind_exc),Trial_num);
% W_I = zeros(length(ind_inh),Trial_num);
% 
% for k=1:Trial_num
% W_E(:,k) = F_synaptic(:,ind_exc)\(I_exc./(0 - F_binary(:,k)));
% W_I(:,k) = F_synaptic(:,ind_inh)\(I_inh./(-75 - F_binary(:,k)));
% end
% 
% hist(a_scale(ind_exc),20)

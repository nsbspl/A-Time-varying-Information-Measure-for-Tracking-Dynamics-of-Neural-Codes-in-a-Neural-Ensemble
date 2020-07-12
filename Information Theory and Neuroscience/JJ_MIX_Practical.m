%clear all
P1 = genpath('R:\MILAD\Results_March_2015\Jitter Sync');%genpath('R:\MILAD\Results_March_2015\J_nCor_ModifiedSigma_40');%genpath('R:\MILAD\Results_March_2015\J_nCor_60percent Sync');%genpath('R:\MILAD\Results_March_2015\J_DynamicClamp');% 

P2 = genpath('R:\MILAD\Results_March_2015\IEEE_CIM');%genpath('C:\Users\Milad Lankarany\Desktop\REcent\IEEE_CIM');
P3 = genpath('C:\Users\Milad Lankarany\Documents\MATLAB\Theoretical');
addpath(P1,P2,P3)
SamplingTime = 0.05;% msec
EndTime = 100;%sec
p = NeuronType(SamplingTime,EndTime,'CD');
tt = 0:p.dt:p.tStop;
% Simulation
midpoint = 15;
a_slow = 60;
a_fast = 85;
% Real data
% midpoint = 70;
% a_slow = 120;%60;%
% a_fast = 170;%85;%

WhichInp = 'mix';
    load Fast_signal, 
    load Slow_signal
    Signal_mix = (midpoint + a_slow*Slow_signal + a_fast*Fast_signal)*1e-3;
    Signal_Inp = Signal_mix;
    
Sigma = (0:5:40)*3;%(0:5:80)*3;%[0;15;60;120;180;240];%0:10:160;%[0;10;50;120];%[10;20;30;40;50;60;70;80;90;100;110];%[0.1; 0.2; 0.3; 0.4; 0.5];
Sigma(1) = 2*3;
%Sigma = [0;0.25;0.5;0.75;1]; % Real data
N = [1;2;5;10;50];%[1;2;5;10];%
CC = [0.05; 0.1; 0.15; 0.2; 0.3; 0.4; 0.5];
tw = 100;
MI = zeros(2,length(Sigma),1);%length(N)
FR = MI;
Err = MI;
MI_std = MI;
CF_mu_M = MI; CF_sem_M = MI;
CF_mu_S = MI; CF_sem_S = MI;
%%
LearningTime = 20/p.dt *1e3;
INP_Learning = Signal_Inp(1:LearningTime);%[Slow_signal(1:LearningTime) Fast_signal(1:LearningTime)];%
Inp_Test = Signal_Inp(LearningTime+1:end); s_orig = (Inp_Test - mean(Inp_Test)); s_orig = s_orig/std(s_orig);
Inp_Test_slow = Slow_signal(LearningTime+1:end); s_orig_slow = (Inp_Test_slow - mean(Inp_Test_slow)); s_orig_slow = s_orig_slow/std(s_orig_slow);
Inp_Test_fast = Fast_signal(LearningTime+1:end); s_orig_fast = (Inp_Test_fast - mean(Inp_Test_fast)); s_orig_fast = s_orig_fast/std(s_orig_fast);
%[f,Power_orig] = PSD(s_orig, p);
for i=3%length(N)%1:length(N)%
    Trial_num = N(i);
    for k = 1:length(Sigma)%1:length(CC)
       %% Learning
        if strcmp(WhichInp,'fast')        
        names = strcat('Ffast',{'_'},'corr',num2str(k),{'_'},'N',num2str(i),{'.mat'}); % this is a cell and the 'save' only accpets char
        C_fast = char(names);
        load(C_fast); 
        F_binary = F_fast;
        elseif strcmp(WhichInp,'slow')       
        names = strcat('Fslow',{'_'},'corr',num2str(k),{'_'},'N',num2str(i),{'.mat'}); % this is a cell and the 'save' only accpets char
        C_slow = char(names);
        load(C_slow); 
        F_binary = F_slow;
        else        
        %names = strcat('Fmix',{'_'},'sigma',num2str(k),{'_'},'N',num2str(i),{'.mat'}); % this is a cell and the 'save' only accpets char
        names = strcat('Fmix',{'_'},'sigma',num2str(k),{'_'},'N',num2str(i),'_jitter03',{'.mat'}); % for jittering case
        C_mix = char(names);
        load(C_mix); 
        F_binary = F_mix;
        end
        
%%
TW = 1;
psth_total = PSTH_(F_binary, p.dt, p.dt);
[F_T, indx_T] = spike_times(F_binary(1:LearningTime,:),p);
psth_Learning = psth_total(1:LearningTime);
FG_mix = KernelPSTH (psth_Learning,TW,p.dt,Trial_num);
TT = find(psth_Learning>0);
% find different classes
theta = FG_mix;
q=2;
Range_ = [min(theta) 0.7*(min(theta) + max(theta)) max(theta)];%linspace(min(theta),max(theta),q+1);%[linspace(min(theta),0.4*(max(theta)+min(theta)),q+1) max(theta)];%[min(theta) linspace(0.4*(max(theta)+min(theta)),max(theta),q+1)];%
indx = [];
for i_=1:q
indx{i_} = find(theta<Range_(i_+1) & theta>=Range_(i_));
%length(indx{i})
end
indx{i_} = sort([indx{i_};find(theta==Range_(i_+1))]);

% Find corresponding spikes
dt = p.dt;
T_sp = [];
for i_=1:q
    aa_ = zeros(size(psth_Learning));
    aa_(indx{i_}) = 1;
    T_sp{i_} = find(psth_Learning.*aa_>0);
    length(T_sp{i_})
end
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
    [sta, stc] = Spike_Triggered_Covariance2(T_sp{i_}*dt, aa(end)*dt, INP_Learning, tw, dt); %INP_Learning
    sta_ = sta - sta(1);
    STA(:,i_) = sta_/(norm(sta_,1)*p.dt);
    
    aa_ = zeros(size(psth_Learning)); 
    aa_(indx{i_}) = 1;
    ps_tt = psth_Learning.*aa_;
end

end
% Find the best coefficients
sig_regul = Signal_mix(1:LearningTime) - midpoint*1e-3;
Inp_M = zeros(length(sig_regul),q);
for i_=1:q
    aa_ = zeros(size(psth_Learning));
    aa_(indx{i_}) = 1;
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
[F_T, indx_T] = spike_times(F_binary(1:LearningTime,:),p);
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
%% Test

psth_Test = psth_total(LearningTime+1:end);
FG_mix = KernelPSTH (psth_Test,TW,p.dt,Trial_num);
TT = find(psth_Test>0);
% find different classes
theta = FG_mix;
q=2;
Range_ = [min(theta) 0.7*(min(theta) + max(theta)) max(theta)];%linspace(min(theta),max(theta),q+1);%[linspace(min(theta),0.4*(max(theta)+min(theta)),q+1) max(theta)];%[min(theta) linspace(0.4*(max(theta)+min(theta)),max(theta),q+1)];%
indx = [];
for i_=1:q
indx{i_} = find(theta<Range_(i_+1) & theta>=Range_(i_));
%length(indx{i})
end
indx{i_} = sort([indx{i_};find(theta==Range_(i_+1))]);
% Find corresponding spikes
dt = p.dt;
T_sp = [];
Inp_Est_M = zeros(length(s_orig),q);
Inp_Est_M = Reconstruct_2((psth_Test-mean(psth_Test)),indx, STA); Inp_Est_M = Inp_Est_M/std(Inp_Est_M);

s_orig_c = s_orig;
Inp_sta = zeros(length(s_orig),1);
Out_Test = psth_total(LearningTime+1:end);
Inp_Est_ = conv(sta,(Out_Test - mean(Out_Test)));%filter(sta,1,(Out_Test - mean(Out_Test)));
Inp_sta = (Inp_Est_(1:length(Out_Test)));% Inp_Est = Inp_Est/std(Inp_Est);
xc = xcorr(s_orig,Inp_sta,200/p.dt);
cnt = 200/p.dt;
[val,ind] = max(xc); 
L2 = length(s_orig) - abs(ind-cnt) + 1;
Inp_sta(1:L2) = Inp_sta(abs(ind-cnt):end); Inp_sta = Inp_sta/std(Inp_sta);
% Out_shuf_ = zeros(size(psth_total)); Out_shuf_(floor(sp_Time_shuf))=1;
% Out_shuf = Out_shuf_(LearningTime+1:end);
% Inp_shuf_ = conv(sta,(Out_shuf - mean(Out_shuf)));
% Inp_shuf = (Inp_shuf_(1:length(Out_Test)))/Trial_num; Inp_shuf = Inp_shuf/std(Inp_shuf);
% Inp_shuf = Inp_shuf(1:length(Inp_Est));
% Fs = 1/(p.dt*1e-3);% Sample time
% L = length(Out_Test);                     % Length of signal
% T = (0:L-1)*(p.dt*1e-3);
% window = 1024;%2^(nextpow2(L)-1);%1256;%
% noverlap = 250;%100;%2000;%2^nextpow2(L)/2;
% nfft = 1024;%2^(nextpow2(L)-1);%1256;%

% [Power_orig,f] = pwelch(Sig_orig/norm(Sig_orig),window,noverlap,nfft,Fs);
% [Power_rec,f] = pwelch(Sig_rec/norm(Sig_rec),window,noverlap,nfft,Fs);

%[f,Power_rec] = PSD(Inp_Est, p);

% figure; 
% semilogx(f,Power_orig,'k')
% hold on,
% semilogx(f,Power_rec,'b')
[Inp_Est_M_slow, Inp_Est_M_fast] = separateSignal(Inp_Est_M, p);
[Inp_sta_slow, Inp_sta_fast] = separateSignal(Inp_sta, p);

mu_MI = [];
std_MI = [];
cf1 = []; cf2 = []; cf3 = []; cf4 = [];
for ll=1:8
    if ll==8
      ty = (ll-1)*(10/p.dt*1e3)+1:length(Inp_Est_M);  
    else
      ty = (ll-1)*(10/p.dt*1e3)+1:ll*(10/p.dt*1e3);
    end
%     [f,Power_noise] = PSD((s_orig(ty) - Inp_Est(ty)), p);
%     [f,Power_shuf] = PSD((s_orig(ty) - Inp_shuf(ty)), p);
%     SNR = Power_orig(ty)./Power_noise(ty);
%     SNR_shuf = Power_orig(ty)./Power_shuf(ty);
%     indx_diff = find(f>50);%find(abs(Power_orig-Power_noise)<=1e-3);
%     SNR(indx_diff(1)+1:end)=0;
%     SNR_shuf(indx_diff(1)+1:end)=0;
%MI_cumulative = cumsum(log2(1+SNR));
%mu_MI = [mu_MI (sum(log2(1+SNR)) - sum(log2(1+SNR_shuf)))];%/sum(Out_Test);
cf1= [cf1 (1 - norm(s_orig_slow(ty) - Inp_Est_M_slow(ty))/norm(s_orig_slow(ty)))];
cf2= [cf2 (1 - norm(s_orig_fast(ty) - Inp_Est_M_fast(ty))/norm(s_orig_fast(ty)))];
cf3= [cf3 (1 - norm(s_orig_slow(ty) - Inp_sta_slow(ty))/norm(s_orig_slow(ty)))];
cf4= [cf4 (1 - norm(s_orig_fast(ty) - Inp_sta_fast(ty))/norm(s_orig_fast(ty)))];
end
% MI_std(1,j,i) = std(mu_MI);
% MI(1,j,i) = mean(mu_MI);
% FR(1,j,i) = sum(Out_Test)/80/Trial_num;
% Err(1,j,i) = norm(s_orig_c - Inp_Est);
CF_mu_M(1,k,i) = mean(cf1); CF_mu_M(2,k,i) = mean(cf2);%1 - norm(s_orig_c - Inp_Est)/norm(s_orig_c);
CF_sem_M(1,k,i) = std(cf1)/sqrt(8); CF_sem_M(2,k,i) = std(cf2)/sqrt(8);
CF_mu_S(1,k,i) = mean(cf3); CF_mu_S(2,k,i) = mean(cf4);%1 - norm(s_orig_c - Inp_Est)/norm(s_orig_c);
CF_sem_S(1,k,i) = std(cf3)/sqrt(8); CF_sem_S(2,k,i) = std(cf4)/sqrt(8);
    end
end
%%

% figure; hold on,
% %for k=1:length(N)
%     errorbar(Sigma(2:end),squeeze(CF_mu_M_Prac(2,2:end,i)),squeeze(CF_sem_M_Prac(2,2:end,i)));
%    % errorbar(Sigma(2:end),squeeze(CF_mu_S(1,2:end,i)),squeeze(CF_sem_S(1,2:end,i)));
figure; hold on,
%for k=1:length(N)
    errorbar(Sigma(2:end)/3,squeeze(CF_mu_M(1,2:end,i)),squeeze(CF_sem_M(1,2:end,i)));
    %errorbar(Sigma(2:end)/3,squeeze(CF_mu_M_Prac(1,2:end,i)),squeeze(CF_sem_M_Prac(1,2:end,i)));
    errorbar(Sigma(2:end)/3,squeeze(CF_mu_S(1,2:end,i)),squeeze(CF_sem_S(1,2:end,i)));
    
figure; hold on,
    errorbar(Sigma(2:end)/3,squeeze(CF_mu_M(2,2:end,i)),squeeze(CF_sem_M(2,2:end,i)));
    %errorbar(Sigma(2:end)/3,squeeze(CF_mu_M_Prac(2,2:end,i)),squeeze(CF_sem_M_Prac(2,2:end,i)));
    errorbar(Sigma(2:end)/3,squeeze(CF_mu_S(2,2:end,i)),squeeze(CF_sem_S(2,2:end,i)));
%end
% save CFmu_M_Prac_corr.mat CF_mu_M_Prac
% save CFsem_M_Prac_corr.mat CF_sem_M_Prac
figure; hold on,
tt_Test = 20e3:p.dt:p.tStop;
plot(tt_Test, s_orig,'k')
plot(tt_Test,Inp_Est_M)
% %plot(tt_Test,Inp_sta,'r')
% 
figure; hold on,
plot(-tw:p.dt:tw,STA)
figure;
plot(-tw:p.dt:tw,sta,'k')

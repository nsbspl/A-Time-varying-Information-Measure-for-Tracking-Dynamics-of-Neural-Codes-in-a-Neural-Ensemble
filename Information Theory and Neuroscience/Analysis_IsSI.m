clear all
P1 = genpath('R:\MILAD\Results_March_2015\Adaptation_SyncEvents');%genpath('R:\MILAD\Results_March_2015\J_nCor_60percent Sync');%genpath('R:\MILAD\Results_March_2015\J_DynamicClamp');% 
addpath(P1)
SamplingTime = 0.05;% msec
EndTime = 100;%sec
p = NeuronType(SamplingTime,EndTime,'CD');
tt = 0:p.dt:p.tStop;
N = [1;2;5;10;50];

percentage_sync = 0.6;
window_sync = 20;% msec
Trial_num = 10;

A_noise = [5;15;25;40]*3;
out_sync = zeros(length(tt),1);
    for k=4%1:length(A_noise)

      
        names = strcat('Fmix',{'_'},'sigma',num2str(A_noise(k)),{'_'},'N',num2str(4),{'.mat'});
        C_mix = char(names);
        load(C_mix); 
        F_binary = F_mix;
dt = p.dt;
q=2;

T_sp = [];
psth_total = PSTH_(F_binary, p.dt, p.dt);
%[F_T, indx_T] = spike_times(F_binary(1:LearningTime,:),p);
[indx_async,indx_sync,sync_event] = find_SyncSP(psth_total,percentage_sync,window_sync,Trial_num,p);
T_sp{1} = indx_async;
T_sp{2} = indx_sync;
out_sync(sync_event) = 1;        
        
        
        %firing_rate(j,k) = sum(F_binary)/EndTime;
        [Pxx_est,fx] = pwelch(out_sync,[],[],[],20e3);%pwelch(F_binary,[],[],[],20e3); % % window ~ length(signal)/10
%         [f1,f2] = find( fx< 2*firing_rate(j,k) & fx> 0.2* firing_rate(j,k) );
%         [p1,p2] = max(Pxx(f1)); peak_psd1(j,k) = 1e8*mean(Pxx(f1(p2))); %1e8*mean(Pxx(f1(p2-3:p2+3)));%
%         [f3,f4] = find(fx<3*firing_rate(j,k) & fx>1.2*firing_rate(j,k));
%         [p3,p4] = max(Pxx(f3)); peak_psd2(j,k) = 1e8*mean(Pxx(f3(p4))); %1e8*mean(Pxx(f3(p4-3:p4+3)));%
%         fr_h1(j,k) = fx(f1(p2)); fr_h2(j,k) = fx(f3(p4));

    end
load Inp_orig    
[Pxx_orig,fx] = pwelch(Inp_orig,[],[],[],20e3);    
indx = find(fx>1&fx<100);
figure;
semilogx(fx(indx),Pxx_orig(indx),'r','LineWidth',3)
hold on,
semilogx(fx(indx),Pxx_est(indx),'k','LineWidth',2)
xlabel ('Freq (Hz)'), ylabel('PSD'), legend('True (input)','Estimated'),title('Using Synchronus Spikes')
[Pxx_2,fx] = pwelch(F_binary(:,5),[],[],[],20e3);
figure;
semilogx(fx(indx),Pxx_orig(indx),'r','LineWidth',3)
hold on,
semilogx(fx(indx),Pxx_2(indx),'b','LineWidth',2)
xlabel ('Freq (Hz)'), ylabel('PSD'), legend('True (input)','Estimated'),title('Using Single Neurons')
%%
th_sync = 0.1;
T_sync = 5; % msec
Trial_num = 50;
TW = 5;
k=2;
A_scale = 20:10:100;
KRT = zeros(length(A_scale),1); SKN = KRT; STD_ = KRT; CV_ = KRT; FR = KRT;
for j=1:length(A_scale)
 

%names = strcat('Fmix',{'_'},'sigmacons',{'_'},'tau',num2str(k),'OU',{'.mat'});
names = strcat('Fmix',{'_'},'scale',num2str(j),{'_'},'tau',num2str(k),'OU',{'.mat'});
        C_mix = char(names);
        load(C_mix); 
        F_binary = F_mix;

psth_total = PSTH_(F_binary, p.dt, p.dt);
[F_T, indx_T] = spike_times(F_binary,p);
psth_Learning = psth_total;
FG_mix = KernelPSTH (psth_Learning,TW,p.dt,Trial_num);
TT = find(psth_Learning>0);
% find different classes
theta = FG_mix(1000:end);
q=2;
Range_ = [min(theta) th_sync*(min(theta) + max(theta)) max(theta)];%linspace(min(theta),max(theta),q+1);%[linspace(min(theta),0.4*(max(theta)+min(theta)),q+1) max(theta)];%[min(theta) linspace(0.4*(max(theta)+min(theta)),max(theta),q+1)];%
indx = [];
for i_=1:q
indx{i_} = find(theta<Range_(i_+1) & theta>=Range_(i_));
%length(indx{i})
end
indx{i_} = sort([indx{i_};find(theta==Range_(i_+1))]);

% Find corresponding spikes
indx_sync = []; amp_sync = [];
[amp_sync,indx_sync] = find_SyncSP(psth_Learning,FG_mix,indx,T_sync,p);
ISI_Sync = [0;diff(indx_sync)]*p.dt;
nbins = 20;
[counts,centers] = hist(ISI_Sync, nbins);

% ISI_Sync = []; ISI_Sync = amp_sync(amp_sync>0);

FR(j) = sum(psth_total)/Trial_num/EndTime;
CV_(j) = std(ISI_Sync)/mean(ISI_Sync);
STD_(j) = std(ISI_Sync);
SKN(j) = skewness(ISI_Sync);
KRT(j) = kurtosis(ISI_Sync);
end
% figure; plot(KRT./FR)
% title('Kurtosis')
% figure; plot(FR)
% title('FR')
% figure; plot(SKN./FR)
% title('Skewness')
% figure; plot(STD_./FR)
% title('std')
% figure; plot(CV_./FR)
% title('Coeff of Var')






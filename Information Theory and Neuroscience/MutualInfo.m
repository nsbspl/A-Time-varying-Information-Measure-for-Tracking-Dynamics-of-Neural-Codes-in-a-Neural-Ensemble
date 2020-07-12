% Calculating the STA
inp = Fast_signal; % change the name accordingly
out = F_mix; % % change the name accordingly
Threshold_sp = -10; % mV
tw = 200; % m-sec The length of sta
spks = spike_binary(out, Threshold_sp); % it gives [0,1]. 1 represents spikes
spk_time = spike_times (spks);
L_W = tw/p.dt;
sta = zeros(2*L_W+1,1);
[sta, stc] = Spike_Triggered_Covariance2(spks, inp, tw, p.dt);
figure; 
plot(-tw:p.dt:tw,sta)
%
cutoff1 = 0.1; cutoff2 = 400;
nfft = 1*1e6;%2*(tw/p.dt);
window = bartlett(2048); noverlap = 1000;
tstep_s = p.dt*1e-3;  %converts to sec
Fs = 1/tstep_s;

stim = inp - mean(inp);
spk = spks - mean(spks);

% Estimated stimulus using different algorithms
IN_sta_ = conv(spk,sta); % KHz
IN_sta_ = IN_sta_(L_W+1:end-L_W);
IN_sta = IN_sta_(1:length(spk)); %IN_sta = IN_sta - mean(IN_sta);

N_sta = stim - IN_sta;

figure; hold on,
plot(stim,'r'), plot(IN_sta,'k')
% Power Spectrums
 [Pstim_True,f] = pwelch(stim,[],[],nfft,Fs); %pwelch(stim,window,noverlap,nfft,Fs);
% [Pstim_sta,f] = pwelch(IN_sta,window,noverlap,nfft,Fs);
% [Pstim_KM,f] = pwelch(IN_KM,window,noverlap,nfft,Fs);
% [Pstim_New,f] = pwelch(IN,window,noverlap,nfft,Fs);
% [Pstim_KK2,f] = pwelch(IN_KK2,window,noverlap,nfft,Fs);
[PN_sta,f] = pwelch(N_sta,[],[],[],Fs); %pwelch(N_sta,window,noverlap,nfft,Fs);
% [PN_KM,f] = pwelch(N_KM,window,noverlap,nfft,Fs);
% [PN_New,f] = pwelch(N_New,window,noverlap,nfft,Fs);
% [PN_KK2,f] = pwelch(N_KK2,window,noverlap,nfft,Fs);
% [Pspk,f] = pwelch(spk,window,noverlap,nfft,Fs);
% [Pstimspk_True, f] = cpsd(stim,spk,window,noverlap,nfft,Fs);
% [Pstimspk_sta, f] = cpsd(IN_sta,stim,window,noverlap,nfft,Fs);
% [Pstimspk_KM, f] = cpsd(IN_KM,stim,window,noverlap,nfft,Fs);
% [Pstimspk_New, f] = cpsd(IN,stim,window,noverlap,nfft,Fs);
% [Pstimspk_KK2, f] = cpsd(IN_KK2,stim,window,noverlap,nfft,Fs);
% Coherence and SNR
disp('computing the coherence and signal-to-noise ratio...');
% Cstimspk_True = (Pstimspk_True).^2./(Pstim_True.*Pspk);
% Cstimspk_sta = (Pstimspk_sta).^2./(Pstim_sta.*Pstim_True);
% Cstimspk_KM = (Pstimspk_KM).^2./(Pstim_KM.*Pstim_True);
% Cstimspk_New = (Pstimspk_New).^2./(Pstim_New.*Pstim_True);
% Cstimspk_KK2 = (Pstimspk_KK2).^2./(Pstim_KK2.*Pstim_True);
% figure; 
% semilogx(f,abs(Cstimspk_sta),'b')
% hold on,
% semilogx(f,abs(Cstimspk_KM),'g')
% semilogx(f,abs(Cstimspk_New),'r')
% semilogx(f,abs(Cstimspk_KK2),'y')
%legend('STA','WK','MSTA','MWK')
%SNRstimspk = 1./(1 - Cstimspk);
% SNR_True = (Cstimspk_True)./(1 - Cstimspk_True);
% SNR_sta = (Cstimspk_sta)./(1 - Cstimspk_sta);%abs((Pstim_sta)./(PN_sta));
% SNR_KM = (Cstimspk_KM)./(1 - Cstimspk_KM);%abs((Pstim_KM)./(PN_KM));
% SNR_New = (Cstimspk_New)./(1 - Cstimspk_New);%abs((Pstim_New)./(PN_New));
% SNR_KK2 = (Cstimspk_KK2)./(1 - Cstimspk_KK2);%abs((Pstim_KK2)./(PN_KK2));
SNR_sta = zeros(size(Pstim_True));
ni=(find(f<=cutoff2 & f>=cutoff1));
SNR_sta_ = Pstim_True./PN_sta; %SNR_sta(SNR_sta<1)=1;
SNR_sta(ni)=SNR_sta_(ni);
% SNR_KM(ni+1:length(f))=0;
% SNR_New(ni+1:length(f))=0;
% SNR_KK2(ni+1:length(f))=0;

figure; 
semilogx(f,log2(1 + abs(SNR_sta)),'b')

% 
% SNRstimspk(ni+1:length(f))=1;
% if ni2>1
%     SNRstimspk(1:ni2-1)=1;
% end
df=f(2)-f(1);
total2=(df*log(SNRstimspk(ni2))/2+df*log(SNRstimspk(ni))/2+df*sum(log(SNRstimspk(ni2+1:ni-1))))/log(2);

1/(4*pi) * df * sum(log2(1+SNR_sta))
%%
% IN_KM_ = conv(spk,h_aa)/(Trial_num); % KHz
% IN_KM_ = IN_KM_(L_W+1:end-L_W);
% IN_KM = IN_KM_(1:length(spk));
% stimest = IN_KM;%fftfilt(h*tstep_s,spk);
%compensates for the delay in the filter
err2_sta = mean( (stim - IN_sta).^2);
err2_KM = mean( (stim - IN_KM).^2);
err2_New = mean( (stim - IN).^2);
err2_KK2 = mean( (stim - IN_KK2).^2);

        %  figure;plot(stimest(nfft+1:length(stimest)));hold on;plot(stim(nfft/2+1+12:length(stim)-nfft/2+12),'r')
        %  pause
          
disp('coding fraction')
cf_sta = 1 - sqrt(err2_sta)/std(stim)
cf_KM = 1 - sqrt(err2_KM)/std(stim)
cf_New = 1 - sqrt(err2_New)/std(stim)
cf_KK2 = 1 - sqrt(err2_KK2)/std(stim)
disp('Mutual Info')
minfo_sta = sum(log2(1 + abs(SNR_sta)))/(p.tStop*1e-3)
minfo_KM = sum(log2(1 + abs(SNR_KM)))/(p.tStop*1e-3)
minfo_New = sum(log2(1 + abs(SNR_New)))/(p.tStop*1e-3)
minfo_KK2 = sum(log2(1 + abs(SNR_KK2)))/(p.tStop*1e-3)
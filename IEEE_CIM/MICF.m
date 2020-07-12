function [cf_rec, minfo_rec] = MICF(spk,stim,stim_rec,p,TotalTime)
%L_W = floor(length(sta_total)/2);
cutoff1 = 0; cutoff2 = 400;
nfft = 5*1e4;%2*(tw/p.dt);
window = bartlett(2048); noverlap = 1000;
tstep_s = p.dt*1e-3;  %converts to sec
Fs = 1/tstep_s;

N_sta = stim - stim_rec;
% Power Spectrums
[Pstim_True,f] = pwelch(stim,window,noverlap,nfft,Fs);
[Pstim_rec,f] = pwelch(stim_rec,window,noverlap,nfft,Fs);
[Pspk,f] = pwelch(spk,window,noverlap,nfft,Fs);
%[PN_sta,f] = pwelch(N_sta,window,noverlap,nfft,Fs);
[Pstimspk_True, f] = cpsd(stim,spk,window,noverlap,nfft,Fs);
[Pstimspk_rec, f] = cpsd(stim_rec,stim,window,noverlap,nfft,Fs);

% Coherence and SNR
Cstimspk_True = (Pstimspk_True).^2./(Pstim_True.*Pspk);
Cstimspk_rec = (Pstimspk_rec).^2./(Pstim_rec.*Pstim_True);
% figure; 
% semilogx(f,abs(Cstimspk_sta),'b')
% hold on,
% semilogx(f,abs(Cstimspk_KM),'g')
% semilogx(f,abs(Cstimspk_New),'r')
% semilogx(f,abs(Cstimspk_KK2),'y')
% legend('STA','WK','MSTA','MWK')
%SNRstimspk = 1./(1 - Cstimspk);
SNR_True = (Cstimspk_True)./(1 - Cstimspk_True);
SNR_rec = (Cstimspk_rec)./(1 - Cstimspk_rec);%abs((Pstim_sta)./(PN_sta));

ni=max(find(f<=cutoff2));
ni2=min(find(f>=cutoff1));
SNR_rec(ni+1:length(f))=0;

% figure; 
% semilogx(f,log2(1 + abs(SNR_sta)),'b')
% hold on,
% semilogx(f,log2(1 + abs(SNR_KM)),'g')
% semilogx(f,log2(1 + abs(SNR_New)),'r')
% semilogx(f,log2(1 + abs(SNR_KK2)),'y')
% legend('STA','WK','MSTA','MWK')


%compensates for the delay in the filter
err2_rec = mean( (stim - stim_rec).^2);


        %  figure;plot(stimest(nfft+1:length(stimest)));hold on;plot(stim(nfft/2+1+12:length(stim)-nfft/2+12),'r')
        %  pause
          
%disp('coding fraction')
cf_rec = 1 - sqrt(err2_rec)/std(stim);
%disp('Mutual Info')
minfo_rec = sum(log2(1 + abs(SNR_rec)))/(TotalTime);

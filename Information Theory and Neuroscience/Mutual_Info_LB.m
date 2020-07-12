function [MI, SNR_sta, f] = Mutual_Info_LB (Inp_true, Inp_est,p)
cutoff1 = 0.1; cutoff2 = 10;
nfft = 1*1e6; Fs = 1/(p.dt*1e-3);
N_sta = Inp_true - Inp_est;
[Pstim_True,f] = pwelch(Inp_true,[],[],[],Fs); %pwelch(stim,window,noverlap,nfft,Fs);
[PN_sta,f] = pwelch(N_sta,[],[],[],Fs); %pwelch(N_sta,window,noverlap,nfft,Fs);


SNR_sta = zeros(size(Pstim_True));
ni=(find(f<=cutoff2 & f>=cutoff1)); 
SNR_sta_ = Pstim_True./PN_sta; %SNR_sta(SNR_sta<1)=1;
SNR_sta(ni)=SNR_sta_(ni);

df=f(2)-f(1);
MI =  df * sum(log2(1+SNR_sta));%1/(4*pi) *
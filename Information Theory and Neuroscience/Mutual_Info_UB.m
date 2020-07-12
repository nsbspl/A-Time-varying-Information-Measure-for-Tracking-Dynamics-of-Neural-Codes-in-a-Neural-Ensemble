function [MI, CRR, f] = Mutual_Info_UB (R_spk,p)
cutoff1 = 0.1; cutoff2 = 10;
nfft = 1*1e6; Fs = 1/(p.dt*1e-3);

[P_ref, f] = cpsd(R_spk(:,1),R_spk(:,2),[],[],[],Fs);
P_ij = zeros(length(P_ref),6);
P_ii = zeros(length(P_ref),3);
k=1;
for i=1:3
    for j=1:3
        if j==i
            [PP, f] = cpsd(R_spk(:,i),R_spk(:,i),[],[],[],Fs);
            P_ii(:,i) = PP;
        else
            [PP, f] = cpsd(R_spk(:,i),R_spk(:,j),[],[],[],Fs);
            P_ij(:,k) = PP;
            k=k+1;
        end
        
    end
    
end
Nom = mean(P_ij,2).^2;
Den = mean(P_ii,2).^2;

CRR = zeros(size(P_ref));
ni=(find(f<=cutoff2 & f>=cutoff1)); 
SNR_sta_ = Nom./Den; %SNR_sta(SNR_sta<1)=1;
CRR(ni)=SNR_sta_(ni);

df=f(2)-f(1);
MI = - df * sum(log2(1-sqrt(CRR)));%1/(4*pi)
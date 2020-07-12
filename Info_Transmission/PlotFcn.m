mil = 0 % layer
k=2; % sigma
i = length(N);
ind = 0*1e3/p.dt+1:1*1e3/p.dt;
if mil==0
        load F_mix_layer_0 % slow signal
        %F_binary = spike_binary(F_mix);
        %load F_mixSig_layer_0 % for the fast signal
        F_binary = F_mix;
    else
        names = strcat('Fmix_2',{'_'},'sigma',num2str(k),{'_'},'layer',num2str(mil),{'.mat'});
        C_mix = char(names);
        load(C_mix); 
        F_binary = F_mix_2;%F_mix_2(:,1:Trial_num);
end
Trial_num = size(F_binary,2);
%% Raster Plot
q=2;
percentage_sync = 0.3;
window_sync = 5;
TW = 1;
psth_total = PSTH_(F_binary, p.dt, p.dt);
psth_ = psth_total(ind);
[indx_async,indx_sync,sync_event,M_S] = find_SyncSP(psth_total,percentage_sync,window_sync,TW,Trial_num,p);
sig_sync = M_S; 
sig_async = ~M_S;
Sig_T{1}=sig_async(ind);
Sig_T{2}=sig_sync(ind);


figure;
%ax(1) = subplot(2,1,1);
col = {'b','r'};%{'k','g','c','y','m','b'};
hold on,
for i=1:Trial_num
    for j=1:q
        jind = (find(F_binary(ind,i).*Sig_T{j}>0)*p.dt)';
        %plot(jind*p.dt,i*ones(length(jind),1),'.','Color',col{j})
        line(repmat(jind,2,1),repmat([i-1;i],1,length(jind)),'color',col{j})
    end
end
hold off
axis([0 length(ind)*p.dt 0 Trial_num])

%figure; plot(ind, Signal_Inp(ind))
% set(gca, 'XTick', [])
% set(gca, 'YTick', [])
% box off
%% PSTHs
fr_neuron = sum(psth_total)/EndTime/Trial_num;
% TW = 5;
% FG_mix = KernelPSTH (psth_total,TW,p.dt,Trial_num);
% scale_ = fr_neuron/mean(FG_mix); % sum(scale_*FG_mix)*p.dt*1e-3/EndTime = sum(psth_total) / Trial_num
% figure; plot(ind,scale_*FG_mix(ind),'k')
% ylim([5 20])

TW = 5; 
FG_mix = KernelPSTH (psth_total,TW,p.dt,Trial_num);
scale_ = fr_neuron/mean(FG_mix);
figure; plot(p.dt*ind,scale_*FG_mix(ind),'k')
%ylim([5 20])
% 
% TW = 100;
% FG_mix = KernelPSTH (F_binary(:,1),TW,p.dt,Trial_num);
% figure; plot(scale_*FG_mix(ind),'k')
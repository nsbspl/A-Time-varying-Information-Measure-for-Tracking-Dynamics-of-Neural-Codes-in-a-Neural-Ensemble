function G = PSTH_(Spike_Binary_, rate, sr)

F_binary = Spike_Binary_; % the spikes of all trials
%rate is in msec
%sr = 0.05; % is always the sampling rate
Trial_num = size(F_binary,2);
rate_ = round(rate/sr);
bin_width_tot = rate;
number_bins_tot = round(length(F_binary)/rate_);
psth_tot = zeros(number_bins_tot,1);
for j=1:number_bins_tot
    if j ==number_bins_tot
        psth_tot(j) = length(find(F_binary((j-1)*rate_+1:end,:)>0));%/(Trial_num*rate*1e-3);%
    else
        psth_tot(j) = length(find(F_binary((j-1)*rate_+1:j*rate_,:)>0));%/(Trial_num*rate*1e-3);%/bin_width_tot;%(j-1)*rate_+1:j*rate_
    end
end
G = psth_tot;
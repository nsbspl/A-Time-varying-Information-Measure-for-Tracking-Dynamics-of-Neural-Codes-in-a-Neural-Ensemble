function sta = Spike_Triggered_Average(sp_binary, stimulus, tw, p)
% sp_binary
%tw = 200; % ms
stimulus = stimulus - mean(stimulus);
length_sta = round(tw/p.dt);
sta = zeros(2*length_sta+1,1);

spike_indx = find(sp_binary>0);
maxindx = spike_indx(end);
indx_int = find(spike_indx>length_sta & spike_indx<maxindx-length_sta);
%sta_resamp = zeros(2*length_stc+1,1);
for j = 1: length(indx_int)
    sta = sta + stimulus(spike_indx(indx_int(j))-length_sta:spike_indx(indx_int(j))+length_sta);
end
sta = sta/length(indx_int);

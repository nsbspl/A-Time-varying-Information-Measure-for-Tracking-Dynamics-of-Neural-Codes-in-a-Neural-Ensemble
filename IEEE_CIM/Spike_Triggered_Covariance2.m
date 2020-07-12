function [sta, stc] = Spike_Triggered_Covariance2(spTime, tEND, stimulus, tw, dt)
% sp_binary
%tw = 200; % ms
rate_=1;
stimulus = stimulus - mean(stimulus);
length_sta = round(tw/dt);
length_stc = length_sta/rate_;
sta = zeros(2*length_sta+1,1);

spike_indx = floor(spTime/dt);%ind(sp_binary>0);
maxindx = floor(tEND/dt);%spike_indx(end);
indx_int = find(spike_indx>length_sta & spike_indx<maxindx-length_sta);
stc = zeros(2*length_stc+1,length(indx_int));
%sta_resamp = zeros(2*length_stc+1,1);
for j = 1: length(indx_int)
    sta = sta + stimulus(spike_indx(indx_int(j))-length_sta:spike_indx(indx_int(j))+length_sta);
    stc_ = stimulus(spike_indx(indx_int(j))-length_sta:spike_indx(indx_int(j))+length_sta);
%     uy = iddata([],stc_,1);
%     ur = resample(uy,length_stc,length_sta);
%     size(ur)
     stc(:,j) = stc_;%ur.u;
end
sta = mean(stc,2);%sta/length(indx_int);
% uy = iddata([],sta,1);
% ur = resample(uy,length_stc,length_sta);

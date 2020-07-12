function sp_shuffle = Shuffle_(Spikes)
F_binary = Spikes;% spike_binary(Voltage);
F_Shuffle = zeros(size(F_binary));
Trial_num = 1;
sp_binary = F_binary;
for i = 1:Trial_num
    %sp_binary (:,i) = F_binary(:,i);
    sp_Total = find(sp_binary(:,i)>0);
    ISI_ = diff(sp_Total);
%hist(ISI_)
%%% shuffling ISI
    num_spike = round(1*length(sp_Total)) - 1;% range of index in integer;
    nn = zeros(num_spike,1);ISI_Shuffle = zeros(size(ISI_));
    nn = floor(num_spike*rand(num_spike,1))+1; % shuffled index of ISI
    for pp = 1:num_spike
        ISI_Shuffle(pp) = ISI_(nn(pp));
    end

    % shuffled spike time
    sptime_shuffle = zeros(num_spike,1);sptime_shuffle(1) = ISI_Shuffle(1);
    for kk = 2:num_spike
        sptime_shuffle(kk) = sptime_shuffle(kk-1) + ISI_Shuffle(kk);
    end
    sptime_shuffle = round(length(F_binary)*(sptime_shuffle/max(sptime_shuffle)));
    %shuffled binary spike trains
    f_shuffle = zeros(length(F_binary),1);
    for kk = 1:length(sptime_shuffle)
       if sptime_shuffle(kk)== 0
         sptime_shuffle(kk) =1;
       end
     f_shuffle(sptime_shuffle(kk))=1;
    end
F_Shuffle(:,i) = f_shuffle;    
end

sp_shuffle = F_Shuffle;
  % plot(tt,f_shuffle)

% [F_T_Shuffle, indx_T] = spike_times(F_Shuffle,p);
% 
% % Raster Plot
% figure;
% hold on,
% %plot(tt,Trial_num*(1/2 + Signal_Total/(2*max(Signal_Total))),'r')
% for i=1:Trial_num
%     plot(F_T_Shuffle(1:indx_T(i),i),i*ones(indx_T(i),1),'*')
% end
% hold off
% xlabel('Time(msec)')
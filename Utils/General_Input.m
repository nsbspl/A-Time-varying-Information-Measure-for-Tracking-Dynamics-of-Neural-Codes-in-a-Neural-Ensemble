
%% Building fast signal
Fast_signal = zeros(L,1);
fast_event = floor(L*(rand(EndTime,1))) + 1; indx_time = [min(fast_event);diff(sort(fast_event))]; 
if min(indx_time)<120/p.dt % to make sure that events are sufficiently far away from each other (> 120 msec)
    indx_time = indx_time + (120/p.dt - min(indx_time));
end
fast_event = cumsum(indx_time);
fast_event(fast_event>L) = L;
Fast_signal(fast_event) = 1;
Fast_signal = conv(Fast_signal,Gz');
Fast_signal = Fast_signal(1:L);

%% Building slow signal
[F_slow] = OUprocess(100,'Gauss',p); 
Slow_signal = F_slow/max(F_slow); Slow_signal = Slow_signal - mean(Slow_signal);
%%
midpoint = 15;
a_slow = 60;
a_fast = 85;
Signal_slow = Slow_signal*a_slow;
Signal_fast = Fast_signal*a_fast;

Eta_noise = OUprocess(5,'Gauss',p); %Eta_noise = Eta_noise/max(Eta_noise);
a_noise = 5;
noise_ = a_noise*Eta_noise;
%%
%p.gNaBar_M = 0;
Signal = (midpoint + Signal_slow + noise_)*1e-3; % Signal_fast 
F_main_1 = singleCompartment_AHP_ML(Signal,p);
F1 = spike_binary(F_main_1(:,1));
sum(F1)/EndTime;
% figure; plot(tt,F_main_1(:,1))
% xlabel('Time (msec)')
% ylabel('Voltage (mV)')
%  figure;
% plot(Signal_slow,'r')
% hold on,
% plot(Signal_fast,'b')
% plot(noise_,'k')

%%
tw = 100;
rate = p.dt;% least time window in msec
sp_Time_L = find(F1>0)*p.dt;
[sta_, stc] = Spike_Triggered_Covariance2(sp_Time_L, sp_Time_L(end), Signal_slow, tw, p.dt);%/Trial_num;
Out_Test = F1;
Inp_Est_ = conv(sta_,(Out_Test - mean(Out_Test)));%filter(sta,1,(Out_Test - mean(Out_Test)));
Inp_Est = (Inp_Est_(1:length(Out_Test))); 
% figure; hold on,
% plot(Signal_slow - mean(Signal_slow),'r')
% plot(Inp_Est,'b')

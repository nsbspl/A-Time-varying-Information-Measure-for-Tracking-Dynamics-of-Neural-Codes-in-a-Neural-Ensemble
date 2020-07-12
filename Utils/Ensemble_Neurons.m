%% Create a population of neurons receiving a common mixed signal + independent noise
Trial_num = 30;
% Signal_Mix = midpoint + Signal_slow + Signal_fast; % pA
Signal_Mix =  zeros(size(Signal_slow)) ; % pA
F_binary = zeros(L,Trial_num);
a_noise = 10; % std of noise (pA)

for k=1:Trial_num  
    Eta_noise = OUprocess(5,'Gauss',p); % std(Eta_noise) = 1 pA
    noise_ = a_noise*Eta_noise;
    Signal = (Signal_Mix + noise_)*1e-3; % converted to nA
    F_main_1 = singleCompartment_AHP_ML(Signal,p);
    F_binary(:,k) = spike_binary(F_main_1(:,1));
end
%% plot
% figure
% for i=1:30
%   hold on
% scatter(tt,F_binary(:,i)*i,'filled');
% 
% end
% title('Spiking Activity of Neurons');
% xlabel('Time (ms)');
% ylabel('Neurons');
% ylim([.5,30.5]);
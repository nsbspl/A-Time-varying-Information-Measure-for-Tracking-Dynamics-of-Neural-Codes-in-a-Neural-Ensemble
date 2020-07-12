psth_total = PSTH_(F_binary, p.dt, p.dt);
aveFR = sum(psth_total)/EndTime/Trial_num;
%%
LearningTime = 7/p.dt *1e3; % in samples
INP_Learning = Signal_Inp(1:LearningTime);%[Slow_signal(1:LearningTime) Fast_signal(1:LearningTime)];%
Inp_Test = Signal_Inp(LearningTime+1:end); 
psth_Learning = psth_total(1:LearningTime);
psth_Test = psth_total(1+LearningTime:end);
%% Finding the Async and Sync time samples
TW = 1;%optw;%
q=2;
dt = p.dt;
T_sp = [];
[indx_async,indx_sync,sync_event] = find_SyncSP(psth_total,percentage_sync,window_sync,TW,Trial_num,p);
T_sp{1} = indx_async;
T_sp{2} = indx_sync;
ps_Total = zeros(length(psth_total),q);
for i = 1:q
    aa_ = zeros(size(psth_total));
    aa_(T_sp{i}) = 1;
    Robs = psth_total.*aa_;
    ps_Total(:,i) = Robs;
end
FG_Ref = zeros(length(psth_total),q);
%% Find the best kernels for each type of spikes
aa_1 = T_sp{1}*p.dt;
aa_2 = T_sp{2}*p.dt;
t = linspace(p.dt,sp_Time(end),1e6); % The total length is 2 M, so the resolution of 1e5 points is equivalent of 1 m-sec which is really good
[yf,tf,optw1,w,c] = sskernel(aa_1);%,t);
[yf,tf,optw2,w,c] = sskernel(aa_2);%,t);
TW_slow = 15;%optw1;%15.07;%optw;
TW_fast = 3;%optw2;%4.33;%;%optw;% optw;
FG_Ref(:,1) = KernelPSTH (ps_Total(:,1),TW_slow,p.dt,Trial_num);
FG_Ref(:,2) = KernelPSTH (ps_Total(:,2),TW_fast,p.dt,Trial_num);
fr1_neuron = sum(ps_Total(:,1))/EndTime/Trial_num; fr2_neuron = sum(ps_Total(:,2))/EndTime/Trial_num;
scale_1 = fr1_neuron/mean(FG_Ref(:,1)); scale_2 = fr2_neuron/mean(FG_Ref(:,2));
%% Learning Process that includes estimating linear and non-linear filters
FG_Total(:,1) = FG_Ref(1:LearningTime,1);%KernelPSTH (Total_PS(:,1),TW_slow,p.dt,Trial_num);
FG_Total(:,2) = FG_Ref(1:LearningTime,2);%KernelPSTH (Total_PS(:,2),TW_fast,p.dt,Trial_num);

filt_sta = stimest2(0,400,FG_Total(:,1),INP_Learning,p.dt,1,2*(tw/p.dt),bartlett(2048),1000);% %(INP_Learning - mean(INP_Learning))
filt_Lin_Async = filt_sta;%(end/2+1:end);
filt_sta = stimest2(0,400,FG_Total(:,2),INP_Learning,p.dt,1,2*(tw/p.dt),bartlett(2048),1000);% %(INP_Learning - mean(INP_Learning))
filt_Lin_Sync = filt_sta;%(end/2+1:end);
L_tw = tw/p.dt;
out_Lin_Async_ = conv(filt_Lin_Async,INP_Learning);%filter(filt_Lin_Async,1,INP_Learning);%
out_Lin_Sync_ = conv(filt_Lin_Sync,INP_Learning);%filter(filt_Lin_Sync,1,INP_Learning);
out_Lin_Async = out_Lin_Async_(L_tw+1:end-L_tw); %out_Lin_Async = out_Lin_Async/std(out_Lin_Async)*std(FG_Total(:,1));
out_Lin_Sync = out_Lin_Sync_(L_tw+1:end-L_tw);

%--- Nonlinearity
num_NL = 20;
NL1 = zeros(num_NL,1); NL2 = NL1;  std_NL1 = NL1; std_NL2 = NL2;
x_lin1 = linspace(min(out_Lin_Async),max(out_Lin_Async),num_NL);
x_lin2 = linspace(min(out_Lin_Sync),max(out_Lin_Sync),num_NL);
for j = 1:num_NL-1
    k1 = find(out_Lin_Async>=x_lin1(j) & out_Lin_Async<x_lin1(j+1));
    NL1(j) = mean(FG_Total(k1,1)); %mean(FG_Total(k1,1));
    std_NL1(j) = std(FG_Total(k1,1));
    k2 = find(out_Lin_Sync>=x_lin2(j) & out_Lin_Sync<x_lin2(j+1));
    NL2(j) = mean(FG_Total(k2,2)); %mean(FG_Total(k2,2));
    std_NL2(j) = std(FG_Total(k2,1));
end
k1 = find(out_Lin_Async==x_lin1(j+1));
NL1(j+1) = mean(FG_Total(k1,1)); %[m,kind] = max(NL1);NL1(kind+1:end)=0;
std_NL1(j+1) = std(FG_Total(k1,1));
k2 = find(out_Lin_Sync==x_lin2(j+1));
NL2(j+1) = mean(FG_Total(k2,2));  NL2(x_lin2<0) = 0; 
std_NL2(j+1) = std(FG_Total(k2,1));

%%%---- Fit the nonlinearity of Async
zeropoint_async = x_lin1(1); k = [];% 
k = find(x_lin1>=zeropoint_async);
xin = x_lin1(k(1:end))';
yout = NL1(k(1:end));
p_async = polyfit(xin,yout,2);
yp_async = nonL_Poly(p_async,x_lin1',zeropoint_async);
figure; hold on,
plot(x_lin1,NL1,'k')
plot(x_lin1,yp_async,'r')
title('Nonlinearity of Asynchronous Spikes')
% Fit the nonlinearity of Synchronous Spikes
xin = x_lin2(1:end-1)';
yout = NL2(1:end-1);
x0 = [150;150;30];
options = optimoptions('lsqcurvefit','TolX',1e-6,'TolFun',1e-6);
p_sync = lsqcurvefit(@nonL_Sync,x0,xin,yout,[],[],options); %polyfit(xin,yout,3);%
% Good values with x = [0.8044; 183.2167, 18.1562]
yp_sync = nonL_Sync(p_sync,x_lin2'); %nonL_Poly(p_sync,x_lin2',0);%
figure; hold on,
plot(x_lin2,NL2,'k')
plot(x_lin2,yp_sync,'r')
title('Nonlinearity of Synchronous Spikes')
% Make thses nonlinerities effects
vq1 = nonL_Poly(p_async,out_Lin_Async,zeropoint_async);%interp1(x_lin1,NL1,out_Lin_Async);
vq2 = nonL_Sync(p_sync,out_Lin_Sync);%interp1(x_lin2,NL2,out_Lin_Sync);
% find the best coef for each class of spikes
ref1 = FG_Ref(1:LearningTime,1);%KernelPSTH (Total_PS(:,1),optw,p.dt,Trial_num);
ref2 = FG_Ref(1:LearningTime,2);%KernelPSTH (Total_PS(:,2),optw,p.dt,Trial_num);
R1 = vq1; coef_Async = (inv(R1'*R1)*R1'*ref1); R2 = vq2; coef_Sync = (inv(R2'*R2)*R2'*ref2);%max(ref2)/max(R2);%
%coef_M = [coef_1;coef_2];%[1;1];% 
% R = [vq1 vq2];
% coef_M = (inv(R'*R)*R'*FG_Test(1:end-1)')
%coef_M = (inv(R'*R)*R'*sum(FG_Ref(1:LearningTime,:),2));  
figure; hold on,
plot(sum(FG_Ref(1:LearningTime,:),2),'k','LineWidth',3) %  sum(FG_Total,2)
plot(vq1*coef_Async + vq2*coef_Sync,'g') %    coef_M
%% Testing Process
% Multiplexed Coding
t_Sync = find(FG_Ref(:,2)>0); % in samples
t_Async = find(FG_Ref(:,2)==0);
ind_Async = t_Async(t_Async>LearningTime)-LearningTime;
ind_Sync = t_Sync(t_Sync>LearningTime)-LearningTime;
FG_Test_Async = FG_Ref(LearningTime+1:end,1);
FG_Test_Sync = FG_Ref(LearningTime+1:end,2);


out_Lin_Async_ = conv(filt_Lin_Async,Inp_Test); out_Lin_Async = out_Lin_Async_(L_tw+1:end-L_tw);
out_Lin_Sync_ = conv(filt_Lin_Sync,Inp_Test); out_Lin_Sync = out_Lin_Sync_(L_tw+1:end-L_tw);
x_lin_test1 = linspace(min(out_Lin_Async),max(out_Lin_Async),num_NL);
x_lin_test2 = linspace(min(out_Lin_Sync),max(out_Lin_Sync),num_NL);
vq1 = nonL_Poly(p_async,out_Lin_Async,zeropoint_async);%interp1(x_lin_test1,NL1,out_Lin_Async);
vq2 = nonL_Sync(p_sync,out_Lin_Sync);%interp1(x_lin_test2,NL2,out_Lin_Sync);
vq_M = vq1*coef_Async + vq2*coef_Sync;%[vq1 vq2]*coef_M;%vq1+vq2;%
vq_M_Total = zeros(length(vq_M),2); %FG_opt_Total = zeros(length(vq_M),2);
vq_M_Total(ind_Async,1) = vq_M(ind_Async); vq_M_Total(ind_Sync,2) = vq_M(ind_Sync);%coef_Sync*vq2(ind_Sync);%

figure; hold on,
plot(FG_Test_Async,'k')
plot(vq_M_Total(:,1),'r')

figure; hold on,
plot(FG_Test_Sync,'k')
plot(vq_M_Total(:,2),'r')

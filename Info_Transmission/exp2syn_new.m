function Gz = exp2syn_new(dt, tEnd, tau_rise, tau_fall)
%tEnd in sec
EndTime = tEnd*1e3;% m-sec
t=(0:dt:EndTime);
% %A=exp(-t/tau1);
% %B=exp(-t/tau2);
% %figure;plot(t,A,'r',t,B,'g');
% tau1 = 0.002;tau2 = 0.02;
% tp = (tau1*tau2)/(tau2-tau1)*log(tau2/tau1);
% factor= -exp(-tp/tau1)+ exp(-tp/tau2);
% f=1/factor;
% Gx= f*(exp(-t/tau2) - exp(-t/tau1));
% 
% tau1 = 0.02;tau2 = 0.20;
% tp = (tau1*tau2)/(tau2-tau1)*log(tau2/tau1);
% factor= -exp(-tp/tau1)+ exp(-tp/tau2);
% f=1/factor;
% Gy=f*(exp(-t/tau2) - exp(-t/tau1));
% 
% 

tau1 = tau_rise;tau2 = tau_fall;
tp = (tau1*tau2)/(tau2-tau1)*log(tau2/tau1);
factor= -exp(-tp/tau1)+ exp(-tp/tau2);
f=1/factor;
Gz=f*(exp(-t/tau2) - exp(-t/tau1));




% figure;plot(t,Gz,'m');
% xlabel('Time (m-sec)');ylabel('Normalized Conductance (nS)');
%legend('tau1=0.002;tau2=0.02','tau1=0.02;tau2=0.2','tau1=0.2;tau2=2.0');
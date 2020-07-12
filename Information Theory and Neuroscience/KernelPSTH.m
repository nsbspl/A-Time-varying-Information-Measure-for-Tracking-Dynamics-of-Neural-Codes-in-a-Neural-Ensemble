function FG = KernelPSTH (psth_tot,TW,res,N)
Tw_new = TW;% in msec   2*1e1*p.dt;% Tw is in ms
bin_width = Tw_new;%/p

t=-5*bin_width:res:5*bin_width;
Gauss_kernel = 1/sqrt(2*pi*(bin_width)^2)*exp(-t.^2/(2*(bin_width)^2)); 
Gauss_kernel = Gauss_kernel/(norm(Gauss_kernel,1)*res);
L_Hkernel = round(length(Gauss_kernel)/2);
Kernel_psth_new = conv(psth_tot,Gauss_kernel);%/(N*res*1e-3); % KHz
Kernel_psth_new = Kernel_psth_new(L_Hkernel:end-L_Hkernel+1);
% uy = iddata([],Kernel_psth_new,1);
% ur = resample(uy,length(psth_tot),length(Kernel_psth_new));
% Kernel_res = ur.u;%/max(ur.u);
%FG = filter(Gauss_kernel,1,psth_tot)/N;%Kernel_psth_new(1:length(psth_tot));
FG = Kernel_psth_new(1:length(psth_tot));






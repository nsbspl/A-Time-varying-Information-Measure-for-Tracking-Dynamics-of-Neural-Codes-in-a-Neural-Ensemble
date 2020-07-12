function [h,Pstimspk] = stimest2(cutoff1,cutoff2,stim,spk,tstep,n,nfft,window,noverlap)

% stimest: function which estimates the stimulus from the spike train
% using the Wiener-Kolmogorov algorithm. It returns the Wiener-Kolmogoro
% filter and the coding fraction which measures the accuracy of the 
% stimulus estimation from the spike train. If the function is called
% without output arguments, the time domain and frequency domain
% characteristics of the estimation are plotted. 
%  
%   [h, tvect, cf] = stimest(stim,spk,tstep)
%
%   where
%       stim = random stimulus
%	spk = spike train
%	tstep = sampling time step (in msec)
%   
%   The function may also be called as follows:
%
%   [h, tvect, cf] = stimest(stim,spk,tstep,n)
%
%   The parameter n determines the subsampling rate (every n-th point is
%   conserved). A reduction in the sampling rate allows  to eliminate 
%   high frequency components outside of the range encoded by the cell and
%   thus eliminates unnecessary noise. 
%
%   Supplementary parameters can be passed by calling:
%
%   [h, tvect, cf] = stimest(stim,spk,tstep,n,nfft,window,noverlap,dflag)
%
%   where the additional parameters control the power spectrum calculation
%   (see matlab psd function) and replace the default values:
%
%       nfft = number of points used for a single fft operation (default: 2048)
%       window = window function (default: bartlett(nfft))
%       noverlap = number of overlapping points per segment (default: 1024)
%       dflag = detrending option flag (default:'none')
%
%   The return parameters are:
%
%   	tvect = vector of time values
%	h = Wiener-Kolmogorov filter
%	cf = coding fraction
%


% if ( (nargin ~= 10) & (nargin ~= 3)  & (nargin ~= 4) )
%   disp(' ');
%   disp('usage1: stimest(stim,spk,tstep) ');
%   disp(' ');
%   disp('usage2: stimest(stim,spk,tstep,n) ');
%   disp(' ');
%   disp('usage3: stimest(stim,spk,tstep,n,nfft,window,noverlap,dflag) ');
%   disp(' ');
%   disp('       for more information type "help psautostim" in the main');
%   disp('       matlab window');
%   disp(' ');
%   return;
% end;


%These parameters are setup with the following simulations in mind:
%a sampling rate of 0.5 msec and a 
%simulation time of 100 sec so that a good estimate of the power
%spectrum is obtained at a resolution of approx. 1Hz. 
% if ( nargin < 7 )
%   nfft = 2048;
%   window = bartlett(nfft);
%   noverlap = 1024;
%   dflag = 'none';
% end;
% 
% if ( nargin == 3 )
%   n = 1;
% end;

stim = stim(:); %converts to column vector if necessary

if ( n > 1) %we need to resample
  disp(' ');
  disp('resampling the data...');
  stim1 = resample(stim,1,n);
  %jb modification to correct problem if original stim length is odd
  stim=stim1(1:length(stim1)-1);
  %stim = stim1;
  clear stim1;
  l_stim = length(stim)
  spk1 = zeros(l_stim,1);
  for k=1:l_stim
    spk1(k,1) = sum(spk((k-1)*n+1:k*n,1));
  end;   
  spk = spk1;
  clear spk1;
  tstep = n*tstep;
end;

%computes the sampling frequency in Hz
tstep_s = tstep*1e-3;  %converts to sec
Fs = 1/tstep_s; %in Hz

%computes and subtracts the mean stimulus value
l_stim = length(stim);
s_stim = sum(stim);
m_stim = s_stim/l_stim;
stim = stim - m_stim;

%computes and subtracts the mean firing rate
spk = spk(:); %converts to column vector if necessary
if length(find(spk==1))>0
spk = spk*Fs; %converts to units of spikes/sec
end
l_spk = length(spk);
s_spk = sum(spk);
m_spk = s_spk/l_spk
spk = spk - m_spk;

disp(' ');
disp('computing the power spectrum of the stimulus...');
[Pstim,f] = pwelch(stim,window,noverlap,nfft,Fs);

disp('computing the power spectrum of the spike train...');
[Pspk,f] = pwelch(spk,window,noverlap,nfft,Fs);

disp('cross-correlating the spike train with the stimulus...');
[Pstimspk, f] = cpsd(stim,spk,window,noverlap,nfft,Fs);


%computes the coherence and signal-to-noise ratio
disp('computing the coherence and signal-to-noise ratio...');
Cstimspk = abs(Pstimspk).^2./(Pstim.*Pspk);
SNRstimspk = 1./(1 - Cstimspk);

%Estimates the transfer function 
disp('computing the Wiener-Kolmogorov filter...');
Tfft_short = Pstimspk./Pspk;

Tfft_long = zeros(nfft,1);
Tfft_long(1:nfft/2+1,1) = Tfft_short(1:nfft/2+1,1);
for k = 2:nfft/2
  Tfft_long(nfft+2-k,1) = conj(Tfft_short(k,1));
end;
%moe's mutual information calculaltion: modified to include lower frequency bound 
ni=max(max(find(f<=cutoff2)));
ni2=min(find(f>=cutoff1));
SNRstimspk(ni+1:length(f))=1;
if ni2>1
    SNRstimspk(1:ni2-1)=1;
end
df=f(2)-f(1);
total2=(df*log(SNRstimspk(ni2))/2+df*log(SNRstimspk(ni))/2+df*sum(log(SNRstimspk(ni2+1:ni-1))))/log(2);
%bialek_mi=(sum(df*log2(SNRstimspk(1:ni))))/2

T = real(ifft(Tfft_long,nfft))*nfft; %Fourier transform with negative phase

%unwraps the filter
tvect = -(nfft/2)*tstep:tstep:(nfft/2)*tstep;
h = zeros(nfft+1,1);
h(1:nfft/2,1) = T(nfft/2+1:nfft,1);
h(nfft/2+1:nfft+1,1) = T(1:nfft/2+1,1);

disp('computing the coding fraction...');
stimest = fftfilt(h*tstep_s,spk);
%compensates for the delay in the filter
err2 = mean( (stimest(nfft+1:length(stimest)) - ...
              stim(nfft/2+1:length(stim)-nfft/2)).^2);

        %  figure;plot(stimest(nfft+1:length(stimest)));hold on;plot(stim(nfft/2+1+12:length(stim)-nfft/2+12),'r')
        %  pause
          
errmax = std(stim);
cf=1 - sqrt(err2)/errmax;
minfo=total2;
%m_info_per_spk=minfo/m_spk
%outf(:,1)=minfo;
%outf(:,2)=m_spk;
%outf(:,3)=m_info_per_spk;
%resutl1=df*mean(Cstimspk(1:20))
%resutl2=df*mean(Cstimspk(40:60))

%figure;plot(f(round(cutoff1)+1:round(cutoff2)),Cstimspk(round(cutoff1)+1:round(cutoff2)),'r');
%  title('Coherence function');
%  xlabel('Frequency [Hz]');
%  ylabel('Coherence value [normalized units]');
%  figure;
% if ( nargout == 0 )
% %looks if the figure 'stimest1' already exists, otherwise creates one
% %and sets it to current
%  % fig_name = 'stimest1';
%  % Figures = get(0,'Chil');
%  % new_fig = 1;
% %  for i=1:length(Figures)
% %    if strcmp(get(Figures(i),'Type'),'figure')
% %      if strcmp(get(Figures(i),'Name'),fig_name)
% %        new_fig = 0;
% %        h_fig = Figures(i);
% %        set(0,'CurrentFigure',h_fig);
% %      end;
% %    end;
% %  end;
% %  if (new_fig == 1)
% %    h_fig = figure('Name',fig_name)
% %    set (gcf,'Position',[1286 33 1260 569]);
% %  end;
% 
% %sets decorations and plots the coherence
% %   subplot(2,2,1);
% %   plot(f(:),Cstimspk(:),'r');
% %   title('Coherence function');
% %   xlabel('Frequency [Hz]');
% %   ylabel('Coherence value [normalized units]');
% %   subplot(2,2,2);
% %   plot(f(:),SNRstimspk(:),'g');
% %   title('Signal-to-Noise Ratio');
% %   xlabel('Frequency [Hz]');
% %   ylabel('Signal-to-Noise Ratio [arbitrary units]'); 
% figure;subplot(2,1,1);
%   plot(tvect,h);
%   title('Filter');
%   xlabel('Time [msec]');
%   ylabel('Filter value [nA]');
%   
%   
%   
%   subplot(2,1,2);
%   
%   plot((0:tstep:10000*tstep),stim(nfft/2+1:nfft/2+1+10000),'k',...
%       (0:tstep:10000*tstep),stimest(nfft+1:nfft+1+10000),'g',...
%      (0:tstep:10000*tstep),spk(nfft/2+1:nfft/2+1+10000)./Fs-4+0.05,'r');
%   title('Stimulus estimation');
%   xlabel('Time [msec]');
%   ylabel('Current [nA]');
% 
%   xlim = get(gca,'XLim');
%   xtot = (xlim(1,2)-xlim(1,1));
%   ylim = get(gca,'YLim');
%   ytot = ylim(1,2) - ylim(1,1);
% 
%   xl1 = [ xlim(1,1)+0.15*xtot xlim(1,1)+0.25*xtot ];
%   yl1 = [ ylim(1,1)+0.1*ytot ylim(1,1)+0.1*ytot ];
%   yl2 = [ ylim(1,1)+0.05*ytot ylim(1,1)+0.05*ytot ];
%   line(xl1,yl1,'Color','k');
%   line(xl1,yl2,'Color','g');
%   sttext = sprintf('stimulus');
%   stesttext = sprintf('est. stimulus');
%   cftext = sprintf('cf = %.2g',cf);
%   text(xl1(1,2)+0.1*xtot,yl1(1,2),sttext);
%   text(xl1(1,2)+0.1*xtot,yl2(1,2),stesttext);
%   text(xlim(1,2)-0.3*xtot,ylim(1,2)-0.1*ytot,cftext);
%   clear h tvect cf;
%   
%   %plot options for ffigure making
%   %bigfig=input('enter 1 for big figure, 0 for not')
%   bigfig=0;
%   if bigfig==1
%       figure
%       set(gcf,'position',[1290 481 1264 458])
%       st_time=input('enter start time in seconds')
%       stp_time=input('enter stop time in sec')
%       st_time=fix((st_time*1e3)/tstep)
%       stp_time=fix((stp_time*1e3)/tstep)
%   
%   plot((st_time*tstep:tstep:stp_time*tstep),stim(nfft/2+1+st_time:nfft/2+1+stp_time),'k',...
%        (st_time*tstep:tstep:stp_time*tstep),stimest(nfft+1+st_time:nfft+1+stp_time),'g',...
%        (st_time*tstep:tstep:stp_time*tstep),spk(nfft/2+1+st_time:nfft/2+1+stp_time)./Fs-4,'r');
% end
% end;

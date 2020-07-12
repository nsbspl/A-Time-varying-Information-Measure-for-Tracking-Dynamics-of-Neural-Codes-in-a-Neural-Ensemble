function bestSigma=FindBestTemporalResolution(f,Visualize)
% find the best temporanl resolution 
% f si spiking activity
% Visualization option for the code

[peakVal,peakIndx]=findpeaks(f);
Thr=max(peakVal);
peakIndx=peakIndx(find(peakVal>=Thr));
peakVal=peakVal(find(peakVal>=Thr));
difrTime(1)=peakIndx(1);

for i = 2:length(peakVal)
    
  difrTime(i)= peakIndx(i)-peakIndx(i-1) ;
    
end
% windowing
Pdfs={};
for i=1:length(peakVal)-1
   windwStrtIndx= peakIndx(i)-round(difrTime(i)/2);
   windwStopIndx= peakIndx(i)+round(difrTime(i+1)/2);
   Pdfs{i} = fitdist(f(windwStrtIndx:windwStopIndx),'Normal');
   sigmas(i)=Pdfs{i}.sigma;
end



if (Visualize)
    maxSigmas=max(sigmas);
     figure
     for i=1:length(peakVal)-1
         w = waitforbuttonpress;
         subplot(211)
         windwStrtIndx= peakIndx(i)-round(difrTime(i)/2);
         windwStopIndx= peakIndx(i)+round(difrTime(i+1)/2);
         
         stem(windwStrtIndx:windwStopIndx,f(windwStrtIndx:windwStopIndx));
         hold on
          scatter(peakIndx(i),f(peakIndx(i)), 'r');
          hold off
         subplot(212)
        
        mu=Pdfs{i}.mu;
        sig=Pdfs{i}.sigma;
        x=mu - 4*maxSigmas:.0001:mu + 4*maxSigmas;
        y = pdf(Pdfs{i},x);
        plot(x,y,'LineWidth',2)
   
end
   end
bestSigma=max(sigmas);
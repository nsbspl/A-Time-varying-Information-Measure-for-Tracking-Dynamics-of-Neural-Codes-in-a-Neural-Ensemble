function [LAE,RMSE,AIC, SSinfo ]= R_LogLikePerformance(yEstimated,yMean, mode,NumFilters,nthist,CellNumb)
 %      Poisson likelihood:      P(s|r) = r^s/s! exp(-r)  
 %   giving log-likelihood:  log P(s|r) =  s log r - r   
%1. for GLM with exponential nonlinearity
[M,N]=size(yEstimated)
SuperSpike=zeros(M,1);
IndSuperSpikw=find(yMean ~= 0);
SuperSpike(IndSuperSpikw)=1;
nsp=sum(SuperSpike);
AIC=zeros(N,1);
SSinfo=zeros(N,1);
RMSE=zeros(N,1);

ratepred_const = nsp/length(yMean);  % mean number of spikes / bin
LL0 = nsp*log(ratepred_const) - length(yMean)*ratepred_const;
if mode==0
    NP=(1+NumFilters);
end
if mode==1
    NP=(1+NumFilters+nthist);
end
if mode==2
    NP=(1+NumFilters+nthist*(CellNumb+1));
end
for i=1:N
     LL_GLM = yMean'*log(yEstimated(:,i)) - sum(yEstimated(:,i));
    LL(i)=LL_GLM;
    SSinfo(i) = (LL_GLM - LL0)/nsp/log(2);
   
    AIC(i) = -2*LL_GLM + 2*NP;
    RMSE(i)=sqrt(mean((yEstimated(:,i)-yMean).^2));
    LAE(i)=(mean(abs(yEstimated(:,i)-yMean)));
end


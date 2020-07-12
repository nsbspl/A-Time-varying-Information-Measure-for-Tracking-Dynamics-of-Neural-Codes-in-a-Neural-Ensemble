function Entropy = entropyFuncTimeIndexWindoL(F_binary,WordLength,pdt)
% This functiono calculate entropy as function of window with length L 
% for each cell
% params:
%   F_binary is spiking activity of the cell
%   WordLength is length of time windodw (structural resolution) to calculate entropy 
%   pdt is the temporal resolution
[SigLen, TotalCells]=size(F_binary); % calculate size of the data 
Entropy=zeros(  SigLen/WordLength,1); % predefinition of variables
 for i = 1 : SigLen  / WordLength
        Wr=F_binary( (i-1)*WordLength+1: (i)*WordLength,:);
        [Pw]=hist(bi2de(Wr'),2^WordLength); % calculate p(w)
         Pw=Pw/length(Wr); % normalize p(w)
       % calculate entropy: E(dt, L)= -1/(dt * L) * sum (p(w) * log(p(w))
        Entropy(i)= -1 * sum(Pw(find (Pw ~= 0)) .* log2(Pw(find (Pw ~= 0))))/(WordLength*pdt)*1000;%
        
end       
         
         
  
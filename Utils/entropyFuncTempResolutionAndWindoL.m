function Entropy = entropyFuncTempResolutionAndWindoL(F_binary, DeltaT,WordLengths,pdt)
% This functiono calculate entropy as function of temporal resolution 
% and window with length L for all cells
% params:
%   F_binary is spiking activity of the cell
%   DeltaT is the temporal resolution
%   WordLength is length of time windodw (structural resolution) to calculate entropy 
%   pdt is the finest temporal resolution
[SigLen, TotalCells]=size(F_binary);
Entropy=zeros( length(DeltaT), length(WordLengths));
Spikes={};
Words={}; 
p=0;
 F_binary=reshape(F_binary,[SigLen * TotalCells,1]);
for worLen = WordLengths
    p=p+1;
    q=0;
for dt = DeltaT
    q=q+1;
    BinLength = (dt/pdt);
    Sp=zeros(length(F_binary) / BinLength,1);
    Wr=zeros(length (Sp)  / worLen,worLen);
   
        for i = 1: length(F_binary)   / BinLength
        if sum( F_binary( (i-1)*BinLength+1: (i)*BinLength)) >=1
            Sp(i) =1;
        end
        end
         for i = 1 : length (Sp )  / worLen
        Wr(i,:)=Sp( (i-1)*worLen+1: (i)*worLen);
         end
         [Pw]=hist(bi2de(Wr),2^worLen); % calculate p(w)
         Pw=Pw/length(Wr); %normalize p(w)
        % calculate entropy: E(Deltat, L)= -1/(Deltat * L) * sum (p(w) * log(p(w))
        %Entropy(q, p)= -1 * sum(Pw(find (Pw ~= 0)) .* log2(Pw(find (Pw ~= 0))))/(worLen*dt)*1000;
        eps=10^-14;
        Entropy(q, p)= -1 * sum(Pw .* log2(Pw+eps))/(worLen*dt)*1000;
end
end
function I=miP(A,B) 
%MI Determines the mutual information of two images or signals
%
%   I=mi(A,B)   Mutual information of A and B, using 256 bins for
%   histograms
%   I=mi(A,B,L) Mutual information of A and B, using L bins for histograms
L=256;

A=double(A); 
B=double(B); 
asyncDist='poisson';

xAsync=linspace( min(A),max(A),L);
dPAsync=fitdist(A,asyncDist );
yAsync = pdf(dPAsync,xAsync);

xSync=linspace( min(B),max(B),L);
dPSync=fitdist(B,asyncDist );
ySync = pdf(dPSync,xSync);

n2 = hist2(A,B,L); 
n2 = n2/sum(n2(:));
I=sum(minf(n2,yAsync'*ySync)); 
% -----------------------
function y=minf(pab,papb)
I=find(papb(:)>1e-12 & pab(:)>1e-12); % function support 
y=pab(I).*log2(pab(I)./papb(I));

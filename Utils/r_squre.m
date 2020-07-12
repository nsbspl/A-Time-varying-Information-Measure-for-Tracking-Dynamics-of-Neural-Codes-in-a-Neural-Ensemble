function r2=r_squre(x,y)
yBar=mean(y);

SSres=sum((x-yBar).^2);
SStot=sum((y-yBar).^2);
r2=1-SSres/SStot;

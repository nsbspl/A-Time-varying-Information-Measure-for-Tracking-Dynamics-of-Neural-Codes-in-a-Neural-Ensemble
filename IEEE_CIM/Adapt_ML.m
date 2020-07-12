function [I] = Adapt_ML(v, G_Adapt, V_Adapt, a)

I = G_Adapt*a*(v-V_Adapt) * 1e-3;  % mA / cm^2 
 
end
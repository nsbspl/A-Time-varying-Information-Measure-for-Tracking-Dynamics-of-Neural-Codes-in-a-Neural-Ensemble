function [I] = Isyn_Exc(v, G_Exc, V_Exc)

I = G_Exc*(v-V_Exc) * 1e-3;  % mA / cm^2 
 
end
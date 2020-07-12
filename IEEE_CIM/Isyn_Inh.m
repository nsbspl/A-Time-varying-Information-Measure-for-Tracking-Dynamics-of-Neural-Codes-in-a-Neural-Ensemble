function [I] = Isyn_Inh(v, G_Inh, V_Inh)

I = G_Inh*(v-V_Inh) * 1e-3;  % mA / cm^2 
 
end
function [I] = Sub_ML(v, G_Sub, V_Sub, z)

I = G_Sub*z*(v-V_Sub) * 1e-3;  % mA / cm^2 
 
end
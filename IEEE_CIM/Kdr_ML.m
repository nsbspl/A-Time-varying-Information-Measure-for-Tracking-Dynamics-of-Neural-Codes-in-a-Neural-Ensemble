function [i, wInf, tau_w] = Kdr_ML(v, w, gK , EK, Beta_w, Gama_w)

wInf = 0.5*(1+tanh((v - Beta_w)/Gama_w));%alpha_n(j)/(alpha_n(j) + beta_n(j));
tau_w = 1./cosh(0.5*(v - Beta_w)/Beta_w);
i = gK*w*(v - EK) * 1e-3;  % mA / cm^2  
end
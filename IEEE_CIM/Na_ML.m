function [i, minf] = Na_ML(v, gNA, ENA, Beta_m, Gama_m)

minf = 0.5*(1+tanh((v - Beta_m)/Gama_m));
i = gNA * minf*(v - ENA) * 1e-3;  % mA / cm^2

end

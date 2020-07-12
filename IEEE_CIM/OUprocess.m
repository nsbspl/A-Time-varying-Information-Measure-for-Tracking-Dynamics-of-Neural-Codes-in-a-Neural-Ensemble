function [ETA,inp_sig] = OUprocess(tau, p)

ETA = zeros(p.tStop/p.dt,1);
inp_sig = ETA;
N_tau = sqrt(2/tau);
k=1;
for t = p.dt:p.dt:p.tStop
    inp_sig(k+1) = randn;
    Einf = tau*N_tau*inp_sig(k)/sqrt(p.dt);
    ETA(k + 1) = Einf + (ETA(k) - Einf)*exp(-p.dt/tau);
    k = k+1;
end
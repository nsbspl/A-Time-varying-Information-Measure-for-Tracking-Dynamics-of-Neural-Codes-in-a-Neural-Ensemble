function [ETA,inp_sig] = OUprocess(tau, stimDist, p)

ETA = zeros(p.tStop/p.dt,1);
inp_sig = ETA;
N_tau = sqrt(2/tau);
k=1;
for t = p.dt:p.dt:p.tStop
    switch stimDist
        case 'exp'
            inp_sig(k+1) = random('exp',1);
        case 'Poisson'
            inp_sig(k+1) = random('Poisson',1);
        case 'Gauss'
            inp_sig(k+1) = randn;
        case 'InvGauss'
            inp_sig(k+1) = random('inversegaussian',1,0.2);
        otherwise
            disp('in make_toy_data: unrecognized probability.')
    end
    
    Einf = tau*N_tau*inp_sig(k)/sqrt(p.dt);
    ETA(k + 1) = Einf + (ETA(k) - Einf)*exp(-p.dt/tau);
    k = k+1;
end
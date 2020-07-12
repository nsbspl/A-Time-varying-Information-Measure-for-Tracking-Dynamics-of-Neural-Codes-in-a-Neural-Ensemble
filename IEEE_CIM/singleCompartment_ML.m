function yy = singleCompartment_ML(Iin_nA, p)
    yy = zeros(p.tStop/p.dt,4);
    v = -70;%p.Vrest;
    n = 0.25*1e-4;%rand;
 
    yy(1,1) = v;
    yy(1,2) = n;
    yy(1,3) = 0;
    yy(1,4) = 0;
    % convert area to cm^2
    A = p.A * 1e-8;    % convert to cm^2
    %Iin_mA = I_inj(t, p) * 1e-6; % convert to mA
    % all currents are mA / cm^2 
    % update IL, the leak current
    k=1;
    for t = p.dt:p.dt:p.tStop
    
    IL = I_Leak(v, p.rL_M, p.EL_M); 
    % update Ikdr an Na
    [iKdr, nInf, taun] = Kdr_ML(v, n, p.gKdrBar_M, p.EK_M, p.V3, p.V4); 
    [iNa,m_inf] = Na_ML(v,p.gNaBar_M, p.ENa_M, p.V1, p.V2);
    dndt = (nInf - n)*p.Phi / taun;
    dvdt = (Iin_nA(k)*1e-6 / A - IL - iKdr - iNa) / (p.cm_M * 1e-3); % dvdt is A / F    % scale_noise*I_noise
    
%     nzn = randn(1);
%     dI_noise = - 1/tau_noise *I_noise + s_n*normn*nzn;
    v = v + p.dt*dvdt;
    n = n + p.dt*dndt;%nInf + (n-nInf)*exp(-p.dt/taun);%;

    %I_noise = I_noise + p.dt*dI_noise;
    yy(k+1,1) = v;
    yy(k+1,2) = n;
    yy(k+1,3) = dvdt;
    yy(k+1,4) = IL;
    k=k+1;
    end
end


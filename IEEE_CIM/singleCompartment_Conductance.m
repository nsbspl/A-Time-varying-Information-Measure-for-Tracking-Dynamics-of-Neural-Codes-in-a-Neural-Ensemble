function yy = singleCompartment_Conductance(I_bias, gExc2, gInh2, Inoise_nA, p)
    yy = zeros(p.tStop/p.dt,5);
    v = -70;%p.Vrest;
    n = 0.25*1e-4;%rand;
    g_exc = rand*1e-0;
    g_inh = rand*1e-0;
    m = 0;
    z = 0; tauZ = p.tauZ;% 40; % m-sec
    alpham=0.005; betam=35; gammam=2;g_M=0;
    alphaz=1/tauZ; betaz=0; gammaz=2; g_AHP=p.g_AHP;% 25;
    exc_avg = p.exc_avg; s_exc = 0;
    tau_exc = 3; % ms
    inh_avg = p.inh_avg; s_inh = 0;
    tau_inh = 10;  I_noise = 0;
    tau_noise = 5; 
    s_n = 15; scale_noise = p.scale_noise;
    normn = sqrt(2/tau_noise);
    yy(1,1) = v;
    yy(1,2) = n;
    yy(1,3) = g_exc;
    yy(1,4) = g_inh;
    yy(1,5) = m;
    yy(1,6) = z;
    % convert area to cm^2
    A = p.A * 1e-8;    % convert to cm^2
    %Iin_mA = I_inj(t, p) * 1e-6; % convert to mA
    % all currents are mA / cm^2 
    % update IL, the leak current
    k=1;
    for t = p.dt:p.dt:p.tStop
    nze = 0.1*randn;
    nzi = 0.1*randn;
    IL = I_Leak(v, p.rL_M, p.EL_M); 
    % update Ikdr an Na
    [iKdr, nInf, taun] = Kdr_ML(v, n, p.gKdrBar_M, p.EK_M, p.V3, p.V4); 
    iNa = Na_ML(v,p.gNaBar_M, p.ENa_M, p.V1, p.V2);
    dndt = (nInf - n)*p.Phi / taun;
    dvdt = (Inoise_nA(k) - 1e-3*( -I_bias + g_exc*(v-0) + g_inh*(v-p.EL_M) + gExc2(k)*(v-0) + gInh2(k)*(v-(-80)) + g_M*m*(v-p.EK_M) + g_AHP*z*(v-p.EK_M)) - IL - iKdr - iNa) * (1e-6 / A) / (p.cm_M * 1e-3); % dvdt is A / F    % scale_noise*I_noise
    dg_exc= -1/tau_exc *(g_exc-exc_avg)+s_exc*nze; %g_inh*(v-p.EL_M)
    dg_inh= -1/tau_inh *(g_inh-inh_avg)+s_inh*nzi;
    dmdt = alpham*(1/(1+exp(-(v+betam)/gammam))-m);
    dzdt = alphaz*(1/(1+exp(-(v+betaz)/gammaz))-z);
%     nzn = randn(1);
%     dI_noise = - 1/tau_noise *I_noise + s_n*normn*nzn;
    v = v + p.dt*dvdt;
    n = n + p.dt*dndt;%nInf + (n-nInf)*exp(-p.dt/taun);%;
    g_exc = 0;%g_exc + p.dt*dg_exc;%
    g_inh = 0;%g_inh + p.dt*dg_inh;% 
    m = m + p.dt*dmdt;
    z = z + p.dt*dzdt;
    %I_noise = I_noise + p.dt*dI_noise;
    yy(k+1,1) = v;
    yy(k+1,2) = n;
    yy(k+1,3) = m;
    yy(k+1,4) = z;
    k=k+1;
    end
end


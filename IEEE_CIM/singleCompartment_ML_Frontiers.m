function yy = singleCompartment_ML_Frontiers(Iin_nA, p)
    yy = zeros(p.tStop/p.dt,4);
    v = -70;%p.Vrest;
    w = 0.25*1e-4;%rand;
    a = 1e-4*rand;
    z = 1e-4*rand;
    g_exc = rand*1e-0;
    g_inh = rand*1e-0;
    exc_avg = p.exc_avg; s_exc = 0;
    tau_exc = 3; % ms
    inh_avg = p.inh_avg; s_inh = 0;
    tau_inh = 10;  I_noise = 0;
    tau_noise = 5; 
    s_n = 15; %scale_noise = p.scale_noise;
    normn = sqrt(2/tau_noise);
    yy(1,1) = v;
    yy(1,2) = w;
    yy(1,3) = g_exc;
    yy(1,4) = g_inh;
    % convert area to cm^2
    A = p.A * 1e-8;    % convert to cm^2
    %Iin_mA = I_inj(t, p) * 1e-6; % convert to mA
    % all currents are mA / cm^2 
    % update IL, the leak current
    k=1;
    for t = p.dt:p.dt:p.tStop
    nze = 0.1*randn;
    nzi = 0.1*randn;
    %-- update Sodium, Potassium and Leak currents
    IL = I_Leak(v, p.rL_M, p.EL_M); 
    [IKdr, wInf, tau_w] = Kdr_ML(v, w, p.gKdrBar_M, p.EK_M, p.Beta_w, p.Gama_w); 
    INa = Na_ML(v,p.gNaBar_M, p.ENa_M, p.Beta_m, p.Gama_m);
    %-- EXC and Inh synaptic conductances
    G_Exc = 0; V_Exc = 0;
    G_Inh = 0; V_Inh = -80;
    I_exc = Isyn_Exc(v, G_Exc, V_Exc);
    I_inh = Isyn_Inh(v, G_Inh, V_Inh);
    %-- Adaptation and Sub currents
    Beta_z = p.Beta_z; Gama_z = p.Gama_z; tau_z = p.tau_z;
    Beta_a = p.Beta_a; Gama_a = p.Gama_a; tau_a = p.tau_a;
    Iadapt = Adapt_ML(v,p.gAdapt_M, p.EK_M, a);
    Isub = Sub_ML(v,p.gSub_M,p.ESub_M, z);
    %-- update the dynamics of the time-varying states (kinematics)
    dzdt = ( 1/( 1+exp((Beta_z-v)/Gama_z) ) - z )/tau_z;
    dadt = ( 1/( 1+exp((Beta_a-v)/Gama_a) ) - a )/tau_a;
    dwdt = (wInf - w)*p.Phi / tau_w;
    dvdt = (Iin_nA(k)*1e-6 / A - I_exc - I_inh - Isub - Iadapt - IL - IKdr - INa) / (p.cm_M * 1e-3); % dvdt is A / F    % scale_noise*I_noise
    %-- update the states 
    v = v + p.dt*dvdt;
    w = w + p.dt*dwdt;%nInf + (n-nInf)*exp(-p.dt/taun);%;
    a = a + p.dt*dadt;
    z = z + p.dt*dzdt;
    %-- save the states
    yy(k+1,1) = v;
    yy(k+1,2) = w;
    yy(k+1,3) = z;
    yy(k+1,4) = a;
    
    k=k+1;
    
    end
end


function Isyn = SynapticFun_Modified (Tspk_1sec,Con_Matrix,E_syn,g_syn,v,Tbl_Exp,t_t,dt)

T_SPK = Tspk_1sec; % The history of the latest 1-sec spikes for all the neurons (active neurons as well as the connected Inactive neurons)

% Tbl_Exp; % the table of synaptic waveform of length 1 sec

t = t_t; % the present time instant

%--- Descriptions of the input: Vol_state,Con_Matrix,E_syn,g_syn are vectors of size N by 1, N is the network size
% Vol_state: includes 0 and 1 representing no-spike and spike of neighboring neurons
% Con_Matrix: includes 0 and 1 representing whether or not the neighboring neurons are connected
% E_syn: the values of sysnaptic reversal potential (depend on the type of synapses)
% g_syn: the strength of synaptic weights
% v: membrane potential of the main neuron at time t
% S_: the history (matrix form) of the synpatic function
% S_generic: is the vector of all generic synaptic function that will be
% added to the S_ if a corresponding neighboring neuron fires
% k: t - tspk, where tspk is the spike time (sample) of all neighboring
% neurons
Isyn = zeros(size(Con_Matrix));
dec_vec = Con_Matrix;% decision vector (or Matrix) --> N*N & The S matrix has a structure similar to dec_vec
% Active synapses are equal to 1.
L = length(Tbl_Exp); N = size(Con_Matrix,2);

[is,js] = find(dec_vec == 1); % (is,js) are the coordinates whenever an spike happen

for i = 1:length(is)
    
    if E_syn(is(i),js(i)) == 0 % exc
          F_Exp = Tbl_Exp(:,1);
      elseif E_syn(is(i),js(i)) == -80 % inh
          F_Exp = Tbl_Exp(:,2);
    end
    ind_true = logical(T_SPK(is(i),:));%find(T_SPK(i,:)>0); % 0s belong to no spike
    indx = (t - T_SPK(is(i),ind_true))/dt;
    indx = 1 + floor(indx);
    if  sum(ind_true) ~= 0
        syn_val = sum(F_Exp(indx));
        Isyn(is(i),js(i)) = syn_val*g_syn(is(i),js(i)).*(E_syn(is(i),js(i)) - v(js(i)));
    else
        Isyn(is(i),js(i)) = 0;
    end
    
end

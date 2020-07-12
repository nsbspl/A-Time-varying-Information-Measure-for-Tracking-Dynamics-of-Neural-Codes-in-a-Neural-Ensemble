function [Isyn,S] = SynapticFun (Vol_state,Con_Matrix,E_syn,g_syn,v,S_,S_generic,k);

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
S = S_;
Isyn = zeros(size(Con_Matrix));
dec_vec = repmat(Vol_state,1,size(Con_Matrix,2)).*Con_Matrix;% decision vector (or Matrix) --> N*N & The S matrix has a structure similar to dec_vec
% Active synapses are equal to 1.
L = length(S_generic); N = size(Con_Matrix,2);
% iL = find(k>L);
% Isyn(iL,:) = 0;
% S(iL,:,:) = zeros(length(iL),N,L);

[is,js] = find(dec_vec == 1); % (is,js) are the coordinates whenever an spike happen

for i = 1:length(is)
    if k(is(i))>L
        S(is(i),js(i),:) = squeeze(S_generic (is(i),js(i),:));
        Isyn(is(i),js(i)) = S(is(i),js(i),1).*(E_syn(is(i),js(i)) - v(js(i))) * g_syn(is(i),js(i));
    else
        S(is(i),js(i),:) = squeeze(S_generic (is(i),js(i),:)) + [squeeze(S_(is(i),js(i),k(is(i))+1:L)); zeros(k(is(i)),1)];
        Isyn(is(i),js(i)) = S(is(i),js(i),1).*(E_syn(is(i),js(i)) - v(js(i))) * g_syn(is(i),js(i));
    end
end

dec_vec_2 = repmat(~Vol_state,1,size(Con_Matrix,2)).*Con_Matrix;% This correspond to connected but no spike
[iz,jz] = find(dec_vec_2 == 1);
for i = 1:length(iz)
    if k(iz(i))>L
        Isyn(iz(i),jz(i)) = 0;
    else
        Isyn(iz(i),jz(i)) = S(iz(i),jz(i),k(iz(i))).*(E_syn(iz(i),jz(i)) - v(jz(i))) * g_syn(iz(i),jz(i));
    end
end

% 
% if k<=L
%     S = repmat(dec_vec,L,1).* S_generic + [S_(k:L,:); zeros(L-k+1,N)];
% else
%     S = repmat(dec_vec,L,1).* S_generic + [S_(k:L,:); zeros(L-k+1,N)];
% end
% Isyn = S(1,:).*(E_syn - v) * g_syn;



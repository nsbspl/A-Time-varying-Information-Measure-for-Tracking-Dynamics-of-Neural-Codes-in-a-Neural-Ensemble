function F = Reconstruct_2_modified(psth_tot,ind, Filters_)
q =size(ind,2);
L_W = floor(length(Filters_)/2);
F_ = zeros(length(Filters_)+length(psth_tot)-1,q);
for i=1:q
    aa_ = zeros(size(psth_tot));
    %ind_flag = isempty(ind{i});
    aa_(ind{i})= psth_tot(ind{i});
    F_(:,i) = conv(aa_-mean(aa_),Filters_(:,i));
end
FF = sum(F_,2);
F = FF(L_W+1:end-L_W);
F = F(1:length(psth_tot));
    
    
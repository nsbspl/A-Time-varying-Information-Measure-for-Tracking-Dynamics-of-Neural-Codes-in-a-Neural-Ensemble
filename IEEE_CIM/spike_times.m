function [F,indd] = spike_times (V,dt)

%global run
L_max = length(find(V==1)); % Total number of spikes in all trials
F = zeros(L_max,size(V,2));
indd = zeros(size(V,2),1);
for i=1:size(V,2)
    a = find(V(:,i)>0);
    F(1:length(a),i) = dt*a;
    indd(i) = length(a);
end


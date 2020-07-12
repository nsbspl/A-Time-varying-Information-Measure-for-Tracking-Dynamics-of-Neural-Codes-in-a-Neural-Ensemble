function G = spike_binary(V)

thresh = -20;
Trial_num = size(V,2);
L = length(V);
F = zeros(size(V,1)+1,size(V,2));

for i = 1:Trial_num
    V1 = [V(:,i);-90]; V2 = [-90; V(:,i)]; 
    F(:,i) = (V1>thresh & V2<thresh);
end
G = F(2:end,:);
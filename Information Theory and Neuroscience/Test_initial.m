std_m = 20;
x = std_m*randn(1e4,1);
%x = poissrnd(std_m,1e4,1);
nbins = 50;
[H, bins] = hist(x,nbins);
figure; hold on,
hist(x,nbins)
plot(bins,H,'r','LineWidth',3)
% the first trick to get rid of negative robabilities (in Hist)
H(H<=0) = eps;
pr_est = H/length(x);

entropy_theory = 1/2 * log2(2*pi*exp(1)*std_m^2)

% the second trick to normalize to the histogram bins
delta_bins = bins(2)-bins(1);
entropy_emprical = -sum(pr_est.*log2(pr_est/delta_bins))

%% Mutual information
% assume that you can estimate the ISI (i.e., x) with an approximation
% which is x_est = x + error
% I call it as p(out/in) because you can get the x_est only ig you are
% given the x

x_est = x + (std_m/10)*randn(size(x));






function RR = Reconstruct(spikes_,resp_,W)

inp_est_ = conv(spikes_,resp_); %
inp_est_ = inp_est_(W+1:end-W); 
inp_est_ = (inp_est_ - mean(inp_est_))/norm((inp_est_ - mean(inp_est_)));
RR = inp_est_;
function i = I_Leak(v, RL, EL)
    
    % determine im
    
gL = 1/RL * 1e3;  % mS / cm^2 
    
    
i = gL * (v - EL) * 1e-3; % mA / cm^2
    

end
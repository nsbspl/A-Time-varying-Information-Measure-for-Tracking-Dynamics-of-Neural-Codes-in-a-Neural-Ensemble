function rmse=RegressionCost(x,EEs,EEa,normSignal_Mix)

wa=x(1);
ws=x(2);
rmse=sqrt(mean(((wa*EEa+ws*EEs)-normSignal_Mix).^2));


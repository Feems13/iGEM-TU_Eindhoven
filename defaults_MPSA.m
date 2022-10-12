function [lb,ub,index_log_parameters] = defaults_MPSA()

index_log_parameters = [12,20,31];

lb = [0.01, 400, 0.006, 10, 1100, 1, 0.05, 0.005, 0.001, 0.0005,...
    0.0005, 0.001, 1090, 0.008, 0.005, 0.020, 0.001, 0.020, 0.001,...
    0.000126, 0.000312, 0.8, 0.4, 0.5, 0.1, 0.2, 0.003, 0.1,...
    0.2, 0.005, 0.000000126, 0.3401, 0, 0.0104, 0.005]; 
ub = [0.01, 400, 0.006, 10, 1100, 1, 0.05, 0.005, 0.001, 0.0005,...
    0.0005, 1000, 1090, 0.008, 0.005, 0.020, 0.001, 0.020, 0.001,...
    1.26, 0.000312, 0.8, 0.4, 0.5, 0.1, 0.2, 0.003, 0.1,...
    0.2, 0.005, 1260, 0.3401, 0, 0.0104, 0.005];

end


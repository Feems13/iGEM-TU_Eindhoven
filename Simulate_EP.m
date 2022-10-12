function [f] = Simulate_EP(p,T,Y0,options)

[T,Y]=ode15s( @(t,y) ODEs_IL10(t,y,p,[]) ,T,Y0,options);

%SEAPex_max = max(Y(:,17));
%SEAPex_t05 = T(round(0.5* find(Y(:,17) == SEAPex_max,1)),1)

%IL10_max = max(IL10ex); % The maximum IL10 value 
%t_36 = find(T == 36); % The index of the time at 36 min
%t_max = t(i_max); % The time at which IL10 has its maximum value
%t_max_h = t_max/3600; % The time at which IL10 has its maximum value in hours
[m48_IL10ex, tmax_IL10ex_detect, tmin_IL10ex_detect] = find_characteristics(T,Y(:,17));

f = [m48_IL10ex, tmax_IL10ex_detect, tmin_IL10ex_detect];
    

end
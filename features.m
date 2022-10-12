function [IL10_1, IL10_2, IL10_5, IL10_10, IL10_15, IL10_20, IL10_25, IL10_30, IL10_36] = features_kd(t, IL10ex)
% Input values are time and IL10 excretion calculated with the ODEs_IL10 file
% Output values are the maximal IL10 excretion with its time in seconds and
% hours. And the half of IL10 excretion with its time in seconds and hours

IL10_1 = IL10ex(3600);
IL10_2 = IL10ex(7200);
IL10_5 = IL10ex(18000);
IL10_10 = IL10ex(36000);
IL10_15 = IL10ex(54000);
IL10_20 = IL10ex(72000);
IL10_25 = IL10ex(90000);
IL10_30 = IL10ex(108000);
IL10_36 = IL10ex(129600);



%IL10_max, t_max, t_max_h, IL10_12, t_12, t_12_h

% IL10_max = max(IL10ex); % The maximum IL10 value 
% i_max = i_max); % The time at which IL10 has its maximum value
% t_max_h = t_max/3600; % The time at which IL10 has its maximum value in hours
% find(IL10ex == IL10_max); % The index of the maximum IL10 value
% t_max = t(
% ex_IL10_12 = IL10_max * (1/2); % Expected IL10 half 
% [val, i_12] = min(abs(IL10ex - ex_IL10_12)); % Gives the absolute difference between the expected value of IL10 and te real IL10 values, the one with the minimal difference will be added
% IL10_12 = IL10ex(i_12); % Half time value of IL10
% t_12 = t(i_12); % Half time value of time
% t_12_h = t_12/3600; % half time value of time in hours

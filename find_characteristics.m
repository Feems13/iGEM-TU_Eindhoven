function [m48,tmax,tmin] = find_characteristics(t,x)

%FUNCTION Calculations/find_characteristics(t, x)
%
%  Description:
%   Find Time of maximum value, relative increase of maximum value, and
%   relative increase at 30 minutes for any x described over t(minutes)
%
%  Inputs:
%   [1] t               - time vector (minutes)
%   [2] x               - any corresponding vector, for which the
%                           characteristics will be determined
%
%  Outputs: 
%   [1] m36              - actual IL-10 concentration
%   [2] tmax             - time of maximal increase (after t(1)) - (double)
%   [3] tmin             - time of minimal increase  (after t(1)) - (double)
% 
%  Contact: Eindhoven University of Technology
% 

%% Maximal measurable concentration IL-10
index_max                   = 0.012;
[val,idx]                   = min(abs(x-index_max));
tmax                        = t(idx);

% Maximal fold change and time
index_min                   = 0.000146;
[val,idx]                   = min(abs(x-index_min));
tmin                        = t(idx);

%% 36 hour IL-10 change
%m36_fold                     = x(find(t==36)) / x(1);
t_48                          = 48*3600; % Switch 36 hours to seconds
m48                          = x(find(t==t_48)); % Find IL-10 value at 36 hours


%% Model defaults
% Call this matlabfile at the starting of a ODE solver to set the defaults
clear all;
close all;

%% Configure time
t_step        = 1;               % [sec]
t_end         = 50*3600;          % [sec] +/- 4 hours
T             = 0:t_step:t_end;  % [sec]

%% Define the labels so that the legend entry for ODE 3,5 are given as labels([3,5])
labels = {  'STAT3c',       'RJ2Lp-STAT3c',  'STAT3cp',      'RJ2Lp-STAT3cp', 'PPX', ...
            'PPX-STAT3cp',  'STAT3cpd',     'STAT3npd',     'STAT3np',      'PPN', ...
            'PPN-STAT3np',  'STAT3n',       'mRNAn-IL-10',   'mRNAc-IL-10',   'IL-10er', ...
            'IL-10g',        'IL-10ex',       'Ligand',       'RJ',           'RJL', ...
            'RJ2L',         'RJ2Lp'};
        
labels_parameters = {'KM1', 'KM2', 'KM3', 'KM4', 'Rtot', 'KMH1', 'KTcSTAT3', 'KTnSTAT3cpd', ...
    'KTcmRNA_IL10', 'KTerIL10', 'KTgIL10', 'Kd1', 'Kd2', 'KB1', 'KB2', 'KB3', 'KB4',...
    'KB5', 'KB6', 'KB7', 'KB8', 'KBu1', 'KBu2', 'KBu4', 'KBu5', 'KBu6', 'KBu7', 'KBu8', 'KBu9',...
    'KBu10', 'KBu11', 'KBu12', 'KDn1', 'KDc1', 'KP1'};
        
units_parameters = {   '[nM/s]',  '[nM]',  '[1/s]',   '[nM]', '[nM]', '[-]', '[1/s]', '[1/s]'...
            '[1/s]',  '[1/s]',  '[1/s]',   '[nM]', '[nM]', '[1/(nM*s)]', '[1/(nM*s)]', '[1/(nM*s)]', '[1/(nM*s)]' ...
            '[1/(nM*s)]',  '[1/(nM*s)]',  '[1/(nM*s)]',   '[1/(nM*s)]', '[1/s]', '[1/s]', '[1/s]', '[1/s]', '[1/s]'  ...
            '[1/s]', '[1/s]', '[1/s]', '[1/s]', '[1/s]', '[1/s]', '[1/s]', '[1/s]', '[1/s]'};

%% Defining initial conditions
% The values bellow are in [nM]
Y0_STAT3c   = 1000;
Y0_3zero    = zeros(1,3);
Y0_PPX      = 50;
Y0_4zero    = zeros(1,4);
Y0_PPN      = 60;
Y0_7zero    = zeros(1,7);
Y0_L        = 0.24;       % The optimum value for the ligand is 1/100 the receptor, found via a Bell curve.
Y0_RJ       = 120;
Y0_3zero;
%     1         2-4        5       6-9       10     11-17      18    19     20-22 
Y0 = [Y0_STAT3c, Y0_3zero, Y0_PPX, Y0_4zero, Y0_PPN, Y0_7zero, Y0_L, Y0_RJ, Y0_3zero];
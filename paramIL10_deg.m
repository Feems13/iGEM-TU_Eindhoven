function [P0] = paramIL10_deg()

% Constants
% Michaelis-Menten
KM1  = 0.01;  % [nM/s]
KM2  = 400;   % [nM]
KM3  = 0.006;  % [1/s] (max 0.5)
KM4  = 10;    % [nM]    (tussen 10 en 65.8)
Rtot = 1100;  % [nM] concentration ribosomes

% Hill exponent 
KMH1 = 1;     % [-] Used for mRNA transcription

% Transport
KTcSTAT3     = 0.05;    % [1/s]
KTnSTAT3cpd  = 0.005;   % [1/s]
KTcmRNA_IL10 = 0.001;   % [1/s]
KTerIL10     = 0.0005;  % [1/s]
KTgIL10      = 0.0005;  % [1/s]

% Binding
Kd1 = 0.741;      % [nM]       (range 1000-0.1 nM) % Kd of ligand/receptor binding
Kd2 = 1090;     % [nM]                           % Kd of the dimerization of the receptor units
KB1 = 0.008;    % [1/(nM*s)]
KB2 = 0.005;    % [1/(nM*s)]
KB3 = 0.02;     % [1/(nM*s)]
KB4 = 0.001;    % [1/(nM*s)]
KB5 = KB3;      % [1/(nM*s)]
KB6 = 0.001;    % [1/(nM*s)]
KB7 = 0.0126;   % [1/(nM*s)] ranges from 0.00001 to 0.0004
KB8 = 0.000312; % [1/(nM*s)]

% Unbinding
KBu1  = 0.8;     % [1/s]
KBu2  = 0.4;     % [1/s]
% KBu3  = nan;   % [1/s]  % KBu3 is unused
KBu4  = 0.5;     % [1/s]
KBu5  = 0.1;     % [1/s]
KBu6  = 0.2;     % [1/s]
KBu7  = 0.003;   % [1/s]
KBu8  = KBu5;    % [1/s]
KBu9  = 0.2;     % [1/s]
KBu10 = 0.005;   % [1/s]
KBu11 = KB7*Kd1; % [1/s] 
KBu12 = KB8*Kd2; % [1/s]

% Degradation 
KDn1 = 0;        % [1/s]
KDc1 = 0.0104;   % [1/s] (Kan lopen van 0.0104 tot 0.78333)
t12_iv = 3600; % [s]
KDe1 = log(2)/t12_iv; %[1/s], 0.0001925

% Phosphorylation  
KP1 = 0.005;     % [1/s]

%% Vector with parameters

P0 = [KM1, KM2, KM3, KM4, Rtot, KMH1, KTcSTAT3, KTnSTAT3cpd, ...
    KTcmRNA_IL10, KTerIL10, KTgIL10, Kd1, Kd2, KB1, KB2, KB3, KB4,...
    KB5, KB6, KB7, KB8, KBu1, KBu2, KBu4, KBu5, KBu6, KBu7, KBu8, KBu9,...
    KBu10, KBu11, KBu12, KDn1, KDc1, KDe1, KP1]; %All parameters for P0 IL_10
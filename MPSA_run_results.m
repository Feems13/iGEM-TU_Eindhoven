%Template for Multi-Parametric Sensitivity Analysis
clear all; close all;
%addpath(genpath('...'))
Defaults_deg; %getting the default settings
%% Initializing cluster
cluster = parcluster('local');

%This code is used when going on a cluster with more power
%if cluster.NumWorkers>4
%n_work = cluster.NumWorkers-1;
%addpath(genpath('....'))
%else
n_work = cluster.NumWorkers;
%end
if isempty(gcp('nocreate'))
pool = parpool(cluster,n_work);
end

%Natal van Riel, TU/e
%load(...) %load model info and struct as needed
par0 = paramIL10_deg(); %Determine the p0. Load the info of the model/initial parameterset
options = odeset('RelTol',1e-9,'nonnegative',1);
[lb,ub,l] = defaults_MPSA();

% Simulate default model for nominal parameter values par0
% and calculate feature from selected output(s)
[f0] = Simulate_EP(par0,T,Y0,options);  %Simulate is a function you have to write yourself

%% Sampling
%p_tweak_loc = [12, 20, 31];
p_tweak_loc = [12,20,31]; %define the indexes of the parameters that are included in the analysis

Np = length(p_tweak_loc);  %all = 35 length(par0) % no. of parameters (rate constants, initial conditions) included in analysis
Nd = 10;       % no. of dummies
n_MPSA = 10000;   % no. of parameter ensembles / no. of Monte Carlo simulations

disp('Monte Carlo simulation information:')
disp(['Number samples: ',num2str(n_MPSA)])
disp(['Parameters tweaked: ',num2str(p_tweak_loc)])

%%

% --- Convert parameters in l to log scale%l=[1,2,3,4,5,7,8,9,11,18,20,22,23,24,25];
lb(l) = log10(lb(l));
ub(l) = log10(ub(l));

% --- Latin hypercub sampling
p_tweak_values = lhsdesign_modified(n_MPSA,lb(p_tweak_loc),ub(p_tweak_loc));
dummy_values = lhsdesign_modified(n_MPSA,zeros(Nd,1),ones(Nd,1));

p_tweak_values_dummies = [p_tweak_values,dummy_values];
% --- Convert log back to normal scale
for li=l
    p_tweak_values(:,p_tweak_loc == li) = 10.^(p_tweak_values(:,p_tweak_loc == li));
end

lb(l) = 10.^(lb(l));
ub(l) = 10.^(ub(l));

% scale = lhsdesign(n_MPSA, Np+Nd);   % random uniform distributed parameter sets (values [0,1])
%                                     % (scale is an n_MPSA x Np+Nd matrix)
%                                     
% % Adjust parameter sets for MC to reference values par0:
% % e.g. 0.1*par0(i) < par_MC(i,j) < 10*par0(i)
% par_MC = scale;
% 
% minp = lb;
% maxp = ub;
% 
% 
% %convert the parameters to log scale for biological parameters 
% minp(index_log_parameters) = log10(lb(index_log_parameters)); 
% maxp(index_log_parameters) = log10(ub(index_log_parameters));
% 
% for i=1:length(p_tweak_loc)
%     MPSA_temp_index = p_tweak_loc(i);
%     par_MC(:,i) = (maxp(p_tweak_loc)- minp(p_tweak_loc)).*scale(:,i) + minp(MPSA_temp_index);
% end
% 
% % --- Convert log back to normal scale
% par_MC(:,find(index_log_parameters == p_tweak_loc)) = 10.^(par_MC(:,find(index_log_parameters == p_tweak_loc)));

%Save initial conditions
save(sprintf('init_MPSA_n%d_p%d_%s.mat',n_MPSA,Np,date),'n_MPSA','Np','Nd','p_tweak_values','p_tweak_values_dummies')

%% Monte Carlo simulations
tic

Ni =4; %4 with high amount of sampling
V = NaN(n_MPSA,length(f0));

for s = 1:Ni
    %parfor i=1:N
    parfor j= 1+n_MPSA/Ni*(s-1):n_MPSA/Ni*s
    %parfor j = 1 : n_MPSA
    
        %partemp=par_MC(j,1:Np); %select the Np model parameters from set j
        partemp = par0;
        partemp(p_tweak_loc) = p_tweak_values(j,:);
        %partemp(MPSA_par_index) = par_MC(j,1:Np);
        % Simulate model for partemp and calculate feature of interest from
        % from selected model outputs
        try
        [f] = Simulate_EP(partemp,T,Y0,options);
        % Calculate sensitivity criterion
        % Sum of squared differences between the perturbed and reference output as example
        V(j,:) = (f0-f).^2;
        
        catch ME
            switch ME.identifier
                case 'MATLAB:UndefinedFunction'
                    warning('Function is undefined.  Assigning a value of NaN.');

                case 'MATLAB:scriptNotAFunction'
                    warning(['Attempting to execute script as function. '...
                        'Running script and assigning output a value of 0.']);
                    notaFunction;
                otherwise
                    rethrow(ME)
            end
        end 
    end
    s
    save(sprintf('results_MPSA_n%d_p%d_%s.mat',n_MPSA,Np,date),'V')
end

toc
%% Closing pool
if ~isempty(gcp('nocreate')) && cluster.NumWorkers>4
delete(gcp('nocreate'))
end

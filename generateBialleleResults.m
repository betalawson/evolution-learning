function generateBialleleResults(regenerate)
% This function generates the results for mass-testing the replicator based
% approaches for the bi-allelic seleciton problem

% Specify the population sizes to run using
pop_sizes = 10.^(2:5);

% Specify how many replicates to use in generating stochastic data to fit
% to, for box plotting
N_rep = 500;

% Specify how many generations to simulate for prediction error checking
sim_gen = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% INPUT HANDLING

% Assume not regenerating data unless specifically requested
if nargin < 1
    regenerate = false;
%elseif regenerate
%    userYN = input("All data will be regenerated. Please confirm Y/N: ",'s');
%    if strcmpi(userYN,'y')
%        regenerate = true;
%    else
%        regenerate = false;
%    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% STRONG SELECTION - s = 0.2 (default for allele_basic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% PROBLEM PREPARATION

% Dataset name
data_prefix = 'DATA_BIALL_';

% Load in the basic problem
load('problems_basic.mat','allele_basic');
problem = allele_basic;

% Specify to use silent mode for mass-testing
ELoptions.nlin_silent = true;


%%% RUN THE PROBLEM ACROSS DIFFERENT METHODS FOR DIFFERENT POPULATION SIZES

% List the problem settings in cell arrays
library_orders = { [0,1], [1], [0,1,2], [1] };
symmetric_flags = { false, false, false, true };

% Set up the filenames where data will be saved/loaded from
data_methodnames = {'orders01','orders1','orders012','symm'};

% Count provided number of methods
N_methods = length(data_methodnames);


%%% GENERATE DATA IF IT IS NOT PRESENT (OR IF ASKED)

% Loop over each method to run and handle each independently
for m = 1:N_methods
    
    % Prepare filename
    this_filename = [data_prefix, data_methodnames{m}, '.mat'];
    
    % Specify options for this method
    ELoptions.library_orders = library_orders{m};
    ELoptions.symmetric_payoff = symmetric_flags{m};
    
    % Generate data for this problem if data not present or requested to
    if ~exist(this_filename,'file') || regenerate
        simResults = runSims( problem, ELoptions, @(library,type) extractBialleleParams(libraryToPayoff(library,2), type), pop_sizes, N_rep, sim_gen );
        save(this_filename, 'simResults', 'ELoptions');
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% WEAKER SELECTION - s = 0.05 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% PROBLEM PREPARATION

% Dataset name
data_prefix = 'DATA_BIALL_WEAK_';

% Load in the basic problem
load('problems_basic.mat','allele_basic');
problem = allele_basic;

% Specify to use silent mode for mass-testing
ELoptions.nlin_silent = true;

% Make selection weaker to test this more challenging case]
s = 0.05;
h = 0.8;
problem.fitness = @(X) [1 + s, 1+h*s; 1+h*s, 1] * X;


%%% RUN THE PROBLEM ACROSS DIFFERENT METHODS FOR DIFFERENT POPULATION SIZES

% List the problem settings in cell arrays
library_orders = { [0,1], [1], [0,1,2], [1] };
symmetric_flags = { false, false, false, true };

% Set up the filenames where data will be saved/loaded from
data_methodnames = {'orders01','orders1','orders012','symm'};

% Count provided number of methods
N_methods = length(data_methodnames);


%%% GENERATE DATA IF IT IS NOT PRESENT (OR IF ASKED)

% Loop over each method to run and handle each independently
for m = 1:N_methods
    
    % Prepare filename
    this_filename = [data_prefix, data_methodnames{m}, '.mat'];
    
    % Specify options for this method
    ELoptions.library_orders = library_orders{m};
    ELoptions.symmetric_payoff = symmetric_flags{m};
    
    % Generate data for this problem if data not present or requested to
    if ~exist(this_filename,'file') || regenerate
        simResults = runSims( problem, ELoptions, @(library,type) extractBialleleParams(libraryToPayoff(library,2), type), pop_sizes, N_rep, sim_gen );
        save(this_filename, 'simResults', 'ELoptions');
    end
    
end
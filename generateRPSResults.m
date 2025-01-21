function generateRPSResults(regenerate)
% This function produces the figures in the paper that demonstrate how
% equation learning performs for identifying rock paper scissors style
% dynamics, which naturally oscillate

% Specify the population sizes to run using
pop_sizes = 10.^(2:5);

% Specify how many replicates to use in generating stochastic data to fit
% to, for box plotting
N_rep = 500;

% Specify how many generations to simulate for prediction error checking
sim_gen = 200;

% Dataset name
data_prefix = 'DATA_RPS_';

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


%%% PROBLEM PREPARATION

% Load in the rock-paper-scissors problems
load('problems_RPS.mat','RPS_balanced','RPS_attract','RPS_repel');

% Create a cell array of cell arrays that contains the problem list
problems_list = {RPS_balanced, RPS_attract, RPS_repel};

% Specify the filename "suffixes" for each problem in the master list
problem_filenames = {'bal','atr','rep'};
problem_names = {'Balanced', 'Attracting', 'Repelling'};

% Specify to use silent mode for mass-testing
ELoptions.nlin_silent = true;


%%% RUN THE PROBLEM ACROSS DIFFERENT METHODS FOR DIFFERENT POPULATION SIZES

% List the problem settings in cell arrays
library_orders = { [0,1], [1], [0,1,2] };

% Provide the additional text added to mark which method was used
method_filenames = {'orders01','orders1','orders012'};

% Count provided number of methods and problems
N_methods = length(method_filenames);
N_problems = length(problems_list);


%%% GENERATE DATA IF IT IS NOT PRESENT (OR IF ASKED)

% Loop over each problem in the problem list
for p = 1:N_problems

    % Loop over each method to run and handle each independently
    for m = 1:N_methods
    
        % Prepare filename
        this_filename = [data_prefix, problem_filenames{p},'_',method_filenames{m}, '.mat'];
    
        % Specify options for this method
        ELoptions.library_orders = library_orders{m};
        ELoptions.symmetric_payoff = false;
        ELoptions.force_positive = true;
    
        % Generate data for this problem if data not present or requested to
        if ~exist(this_filename,'file') || regenerate
            simResults = runSims( problems_list{p}, ELoptions, @(library,type) extractRPSParams(libraryToPayoff(library,3), type), pop_sizes, N_rep, sim_gen );
            save(this_filename, 'simResults', 'ELoptions');
        end
   
    end
end
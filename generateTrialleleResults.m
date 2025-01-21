function generateTrialleleResults(regenerate)
% This function generates the results for mass-testing the replicator based
% approaches for the tri-allelic seleciton problem

% Specify the population sizes to run using
pop_sizes = 10.^(2:5);

% Specify how many replicates to use in generating stochastic data to fit
% to, for box plotting
N_rep = 500;

% Specify how many generations to simulate for prediction error checking
sim_gen = 50;

% Dataset name
data_prefix = 'DATA_TRIALL_';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% INPUT HANDLING

% Assume not regenerating data unless specifically requested
if nargin < 1
    regenerate = false;
%elseif regenerate
    %userYN = input("All data will be regenerated. Please confirm Y/N: ",'s');
    %if strcmpi(userYN,'y')
    %    regenerate = true;
    %else
    %    regenerate = false;
    %end
end


%%% PROBLEM PREPARATION

% Load in the tri-allelic selection problems
load('problems_inheritance.mat', 'allele_standard', 'allele_transient', 'allele_persistence');

% Create a cell array of cell arrays that contains the problem list
problems_list = { allele_standard, allele_transient, allele_persistence };

% Specify the filename "suffixes" for each problem in the master list
problem_filenames = {'std', 'tran', 'pers'};

% Specify to use silent mode for mass-testing
ELoptions.nlin_silent = true;


%%% RUN THE PROBLEM ACROSS DIFFERENT METHODS FOR DIFFERENT POPULATION SIZES

% List the problem settings in cell arrays
library_orders = { [0,1], [1], [0,1,2], [1] };
symmetric_flags = { false, false, false, true };

% Set up the filenames where data will be saved/loaded from
method_filenames = {'orders01','orders1','orders012','symm'};

% Count provided number of methods and number of problems
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
        ELoptions.symmetric_payoff = symmetric_flags{m};
    
        % Generate data for this problem if data not present or requested to
        if ~exist(this_filename,'file') || regenerate
            simResults = runSims( problems_list{p}, ELoptions, @(library,type) extractTrialleleParams(libraryToPayoff(library,3), type), pop_sizes, N_rep, sim_gen );
            save(this_filename, 'simResults', 'ELoptions');
        end
   
    end
end


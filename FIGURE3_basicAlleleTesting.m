function FIGURE3_basicAlleleTesting(regenerate)
% This function produces the figures in the paper that demonstrate how
% equation learning performs for identifying selection strength and
% dominance for a bi-allelic inheritance scenario

% Specify the population sizes to run using
pop_sizes = 10.^(2:5);

% Specify whether to save figures or just leave them all displayed
save_figs = false;

% Specify the population level to compare methods for
pop_level = 2;

% Specify how many replicates to use in generating stochastic data to fit
% to, for box plotting
N_rep = 500;

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

% Load in the basic problem - type I and type II selection versions.
% Settings for these will be adjusted as needed
load('problems_basic.mat','allele_basic');
problem = allele_basic;

% Specify to use silent mode for mass-testing
ELoptions.nlin_silent = true;


%%% RUN THE PROBLEM ACROSS DIFFERENT METHODS FOR DIFFERENT POPULATION SIZES

% List the problem settings in cell arrays
library_orders = { [0,1], [1], [0,1,2], [1] };
symmetric_flags = { false, false, false, true };

% Set up the filenames where data will be saved/loaded from
data_prefix = 'DATA_BIALL_';
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
        simResults = runSims( problem, ELoptions, @(K,Ftexts,type) extractBialleleParams(fitnessToPayoff(K,Ftexts,2), type), pop_sizes, N_rep );
        save(this_filename, 'simResults', 'ELoptions');
    end
    
end

% Safety pause
pause(1);


%%% PLOT FIGURE 3A: TRAJECTORY FIT QUALITY

% Load in the data for each method
for m = 1:N_methods
    
    % Prepare filename
    this_filename = [data_prefix, data_methodnames{m}, '.mat'];
    % Load data
    load(this_filename, 'simResults');
       
    % Use the first dataset to prepare storage and check population size
    if m == 1
        
        % Count number of replicates and warn user if different to above
        N_reps_here = size(simResults.l2_mat,3);
        if N_rep ~= N_reps_here
            warning('Number of replicates in data does not match number of replicates specified in data generator. Check data and re-generate if required.');
            N_rep = N_reps_here;
        end
        
        % Prepare results storage
        l2_mat = zeros(N_methods, 4, N_rep);
        
        % Note population size
        N_pop = simResults.pop_sizes(pop_level);
        
    end
    
    % Grab out the data for the requested population level and store
    l2_mat(m, :, :) = simResults.l2_mat(pop_level,:,:);
    
end

% Prepare all text for figure
x_txts = {'p \leq 1', 'p = 1', 'p \leq 2', 'Symm'};
title_txt = ['N = ',num2str(N_pop)];

% Plot the error boxplots
plotRMSEBoxplot( struct('l2_mat', l2_mat), x_txts, title_txt, 'FIG3a_methodRMSE', save_figs, false);



%%% PLOT FIGURE 3B: RECOVERY OF PARAMETERS FOR DIFFERENT POPULATION SIZES

% Load in the data for the symmetric case
load([data_prefix,'symm.mat'], 'simResults');

% Prepare the x tick labels
for n = 1:length(simResults.pop_sizes)
    x_txts{n} = num2str(simResults.pop_sizes(n));
end
x_txts{end+1} = '\infty';

% Get the true values of the parameters and add these to the results struct
A = problem.fitness(eye(2));
simResults.params_true = extractBialleleParams(A,1);

plotParamsBoxplot( simResults, x_txts, '', 'FIG3b_popsParams', save_figs, true );